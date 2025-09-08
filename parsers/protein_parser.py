"""HMDB protein XML parser for extracting association data."""

import uuid
from typing import Dict, Iterator, List, Tuple

from lxml import etree as ET

from utils.cache_helper import load_pickle
from utils.parser_helper import IDHierarchy, XMLParseHelper


class HMDBProteinParser(XMLParseHelper):
    """Parser for HMDB protein XML data."""

    def __init__(self, input_xml: str):
        super().__init__(input_xml)

        self.uniprot2entrezgene = load_pickle("hmdbp_uniprot2entrezgene.pkl") or {}
        self.protein_functions = self._flatten_protein_functions(
            load_pickle("hmdbp_uniprot_protein_functions.pkl") or []
        )
        self.gene_summaries = self._flatten_gene_summaries(
            load_pickle("hmdbp_entrezgene_summaries.pkl") or []
        )
        self.cached_pathway_descriptions = load_pickle("smpdb_pathway_descriptions.pkl") or {}
        self.cached_go_descriptions = load_pickle("hmdbp_go_definitions.pkl") or {}

    def _flatten_protein_functions(self, protein_function_list: List[Dict]) -> Dict:
        """Flatten list of protein function dictionaries."""
        flattened = {}
        for func_dict in protein_function_list:
            if func_dict:
                flattened.update(func_dict)
        return flattened

    def _flatten_gene_summaries(self, gene_summary_list: List[Dict]) -> Dict:
        """Flatten list of gene summary dictionaries."""
        flattened = {}
        for summary_dict in gene_summary_list:
            if summary_dict:
                flattened.update(summary_dict)
        return flattened

    def get_list_of_region_list(self, elem, tag: str) -> List[List[int]]:
        """Extract list of numeric ranges from XML elements."""
        ranges = []
        if elem is None:
            return ranges

        for e in elem.findall(f"hmdb:{tag}", self.namespace):
            if not e.text or "-" not in e.text:
                continue

            parts = e.text.split("-", 1)
            if len(parts) == 2:
                start, end = (part.strip() for part in parts)
                if start.isdigit() and end.isdigit():
                    ranges.append([int(start), int(end)])

        return ranges

    def get_protein_properties(self, protein) -> Dict:
        """Extract protein properties including sequences and domains."""
        props = protein.find("hmdb:protein_properties", self.namespace)
        if props is None:
            return {}

        # basic properties
        residue_num = self.get_text(props, "residue_number")
        molecular_weight = self.get_text(props, "molecular_weight")
        pi = self.get_text(props, "theoretical_pi")
        prot_seq = self.get_text(props, "polypeptide_sequence")

        # gene regions
        tm_elem = props.find("hmdb:transmembrane_regions", self.namespace)
        tm_regions = self.get_list_of_region_list(tm_elem, "region")

        sig_elem = props.find("hmdb:signal_regions", self.namespace)
        sig_regions = self.get_list_of_region_list(sig_elem, "region")

        # Pfam domains
        pfams_elem = props.find("hmdb:pfams", self.namespace)
        pfam_list = []
        if pfams_elem is not None:
            for pfam in pfams_elem.findall("hmdb:pfam", self.namespace):
                pfam_id = self.get_text(pfam, "pfam_id")
                name = self.get_text(pfam, "name")
                if pfam_id:
                    pfam_list.append(
                        {
                            "id": f"PFAM:{pfam_id}",
                            "name": name.lower() if name else None,
                        }
                    )

        return self.remove_empty_none_values(
            {
                "residue_num": int(residue_num) if residue_num else None,
                "molecular_weight": float(molecular_weight) if molecular_weight else None,
                "theoretical_pi": float(pi) if pi else None,
                "transmembrane_region": tm_regions if tm_regions else None,
                "signal_region": sig_regions if sig_regions else None,
                "protein_seq": prot_seq,
                "pfam": pfam_list if pfam_list else None,
            }
        )

    def get_gene_properties(self, protein) -> Dict:
        """Extract gene properties including chromosomal location."""
        props = protein.find("hmdb:gene_properties", self.namespace)
        if props is None:
            return {}

        chrom_loc_raw = self.get_text(props, "chromosome_location")
        locus = self.get_text(props, "locus")
        gene_seq = self.get_text(props, "gene_sequence")

        # parse chromosomal location
        chrom_loc = None
        if chrom_loc_raw:
            candidate = chrom_loc_raw.split(":", 1)[-1] if ":" in chrom_loc_raw else chrom_loc_raw
            try:
                chrom_loc = int(candidate)
            except ValueError:
                chrom_loc = candidate.strip()

        return self.remove_empty_none_values(
            {
                "chromosomal_location": chrom_loc,
                "locus": locus,
                "gene_sequence": gene_seq,
            }
        )

    def get_publications(self, protein) -> Dict:
        """Extract publication references from protein."""
        pmids = []
        ref_elem = protein.find("hmdb:general_references", self.namespace)
        if ref_elem is not None:
            for ref in ref_elem.findall("hmdb:reference", self.namespace):
                pmid_elem = ref.find("hmdb:pubmed_id", self.namespace)
                if pmid_elem is not None and pmid_elem.text:
                    try:
                        pmids.append(f"PMID:{int(pmid_elem.text.strip())}")
                    except ValueError:
                        continue

        if pmids:
            pmids = sorted(set(pmids))
            return {"pmid": pmids, "category": "biolink:Publication"}

        return {}

    def lookup_uniprot_data(self, uniprot_id: str) -> Tuple[str, str]:
        """Lookup UniProt description and Entrez gene mapping."""
        description = self.protein_functions.get(uniprot_id, {}).get("description")
        entrez_map = self.uniprot2entrezgene.get(uniprot_id, {})
        entrez_curie = entrez_map.get("gene_id")
        return description, entrez_curie

    def lookup_entrez_data(self, entrez_curie: str) -> Tuple[str, str] | Tuple[None, None]:
        """Lookup Entrez gene ID and description."""
        if not entrez_curie or ":" not in entrez_curie:
            return None, None

        entrez_id = entrez_curie.split(":", 1)[1]
        gene_description = self.gene_summaries.get(entrez_id, {}).get("description")
        return entrez_id, gene_description

    def build_protein_node(self, protein) -> Dict:
        """Build standardized protein node with all available information."""
        name = self.get_text(protein, "gene_name")
        full_name = self.get_text(protein, "name")

        # synonyms
        synonyms_elem = protein.find("hmdb:synonyms", self.namespace)
        synonyms = self.get_list(synonyms_elem, "synonym") if synonyms_elem is not None else []

        # primary ID and xrefs
        primary_id, xrefs = IDHierarchy.get_protein_primary_id(protein, self)

        # protein function information
        prot_func = self.get_text(protein, "general_function")
        prot_spec_func = self.get_text(protein, "specific_function")
        prot_type = self.get_text(protein, "protein_type")

        # cellular components
        cellular_com_elem = protein.find("hmdb:subcellular_locations", self.namespace)
        cellular_components = (
            self.get_list(cellular_com_elem, "subcellular_location")
            if cellular_com_elem is not None
            else []
        )

        # PDB structures
        pdb_elem = protein.find("hmdb:pdb_ids", self.namespace)
        pdbs = self.get_list(pdb_elem, "pdb_id") if pdb_elem is not None else []
        if pdbs:
            xrefs["pdb"] = [f"PDB:{pdb}" for pdb in pdbs]

        # protein and gene properties
        prot_props = self.get_protein_properties(protein)
        gene_props = self.get_gene_properties(protein)

        # add Pfam domains to xrefs
        if prot_props.get("pfam"):
            xrefs["pfam"] = prot_props["pfam"]

        # lookup external data
        uniprot_id = self.get_text(protein, "uniprot_id")
        prot_description, entrez_curie, gene_description = None, None, None

        if uniprot_id:
            prot_description, entrez_curie = self.lookup_uniprot_data(uniprot_id)
            if entrez_curie:
                entrez_id, gene_description = self.lookup_entrez_data(entrez_curie)
                if entrez_id:
                    xrefs["entrezgene"] = f"NCBIGene:{entrez_id}"

        # complete protein node
        protein_node = {
            "id": primary_id,
            "name": name,
            "full_name": full_name.lower() if full_name else None,
            "synonyms": synonyms,
            "description": prot_description,
            "function": prot_func,
            "specific_function": prot_spec_func,
            "residue_num": prot_props.get("residue_num"),
            "molecular_weight": prot_props.get("molecular_weight"),
            "theoretical_pi": prot_props.get("theoretical_pi"),
            "transmembrane_region": prot_props.get("transmembrane_region"),
            "signal_region": prot_props.get("signal_region"),
            "protein_seq": prot_props.get("protein_seq"),
            "chromosomal_location": gene_props.get("chromosomal_location"),
            "locus": gene_props.get("locus"),
            "gene_seq": gene_props.get("gene_sequence"),
            "gene_description": gene_description,
            "category": "biolink:Protein",
            "protein_type": prot_type.lower() if prot_type else None,
            "cellular_component": cellular_components,
            "xrefs": xrefs,
        }

        return self.remove_empty_none_values(protein_node)

    def parse_protein_pathway(self) -> Iterator[Dict]:
        """Parse protein-pathway associations."""
        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}protein':
                protein = elem

                # protein node
                subject_node = self.build_protein_node(protein)
                if not subject_node:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                # publication references
                publication = self.get_publications(protein)

                # pathways
                pathways_elem = protein.find("hmdb:pathways", self.namespace)
                if pathways_elem is None:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                for pathway in pathways_elem.findall("hmdb:pathway", self.namespace):
                    pw_name = self.get_text(pathway, "name")
                    if not pw_name:
                        continue

                    kegg_map = self.get_text(pathway, "kegg_map_id")

                    # SMPDB information
                    cache_entry = self.cached_pathway_descriptions.get(pw_name.lower(), {})
                    smpdb_id = cache_entry.get("smpdb")
                    description = cache_entry.get("description")

                    # pathway node
                    object_node = {
                        "id": smpdb_id or (f"KEGG:{kegg_map}" if kegg_map else None),
                        "name": pw_name.lower(),
                        "description": description,
                        "category": "biolink:Pathway",
                        "xrefs": {
                            "smpdb": smpdb_id,
                            "kegg": f"KEGG:{kegg_map}" if kegg_map else None,
                        },
                    }
                    object_node = self.remove_empty_none_values(object_node)

                    # association node
                    association_node = {
                        "category": "biolink:GeneToPathwayAssociation",
                        "predicate": "biolink:participates_in",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "has_evidence": "ECO:0000305",  # manual assertion
                        "agent_type": "biolink:manual_agent",
                        "publication": publication if publication else None,
                    }
                    association_node = self.remove_empty_none_values(association_node)

                    # association ID
                    _id = (
                        f"{subject_node['id'].split(':')[1]}"
                        f"_{association_node['predicate'].split(':')[1]}"
                        f"_{object_node['id'].split(':')[1]}"
                        if "id" in object_node and "id" in subject_node
                        else f"uuid:{str(uuid.uuid4())}"
                    )

                    yield {
                        "_id": _id,
                        "association": association_node,
                        "object": object_node,
                        "subject": subject_node,
                    }

                elem.clear()
                for ancestor in elem.xpath("ancestor-or-self::*"):
                    ancestor.clear()

    def parse_protein_biological_process(self) -> Iterator[Dict]:
        """Parse protein-biological process associations."""
        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}protein':
                protein = elem

                subject_node = self.build_protein_node(protein)
                if not subject_node:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                # publication references
                publication = self.get_publications(protein)

                # GO biological processes
                go_elem = protein.find("hmdb:go_classifications", self.namespace)
                if go_elem is None:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                for go_class in go_elem.findall("hmdb:go_class", self.namespace):
                    category = self.get_text(go_class, "category")
                    if category != "Biological process":
                        continue

                    go_name = self.get_text(go_class, "description")
                    go_id = self.get_text(go_class, "go_id")

                    # GO description from cache
                    go_description = None
                    if go_id and go_id in self.cached_go_descriptions:
                        go_description = self.cached_go_descriptions[go_id].get("description")

                    # biological process node
                    object_node = {
                        "id": go_id,
                        "name": go_name.lower() if go_name else None,
                        "description": go_description,
                        "category": "biolink:BiologicalProcess",
                        "xrefs": {"go": go_id} if go_id else {},
                    }
                    object_node = self.remove_empty_none_values(object_node)

                    # association node
                    association_node = {
                        "category": "biolink:MacromolecularMachineToBiologicalProcessAssociation",
                        "predicate": "biolink:participates_in",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "has_evidence": "ECO:0000305",  # manual assertion
                        "agent_type": "biolink:manual_agent",
                        "publication": publication if publication else None,
                    }
                    association_node = self.remove_empty_none_values(association_node)

                    # association ID
                    _id = (
                        f"{subject_node['id'].split(':', 1)[1]}"
                        f"_{association_node['predicate'].split(':')[1]}"
                        f"_{object_node['id'].split(':', 1)[1]}"
                        if "id" in object_node and "id" in subject_node
                        else f"uuid:{str(uuid.uuid4())}"
                    )

                    yield {
                        "_id": _id,
                        "association": association_node,
                        "object": object_node,
                        "subject": subject_node,
                    }

                elem.clear()
                for ancestor in elem.xpath("ancestor-or-self::*"):
                    ancestor.clear()
