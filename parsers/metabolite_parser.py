"""HMDB metabolite XML parser for extracting association data."""

import uuid
from typing import Dict, Iterator, List

from lxml import etree as ET

from utils.cache_helper import load_pickle
from utils.parser_helper import (
    IDHierarchy,
    MolecularWeightExtractor,
    OrganismClassifier,
    PropertyExtractor,
    ReferenceExtractor,
    XMLParseHelper,
)


class HMDBMetaboliteParser(XMLParseHelper):
    """Parser for HMDB metabolite XML data."""

    def __init__(self, input_xml: str):
        super().__init__(input_xml)
        self.reference_extractor = ReferenceExtractor()

        self.cached_taxon_info = load_pickle("original_taxon_name2taxid.pkl") or {}
        self.cached_disease_info = load_pickle("original_disease_name2id.pkl") or {}
        self.cached_protein_functions = self._flatten_protein_functions(
            load_pickle("uniprot_protein_functions.pkl") or []
        )
        self.cached_pathway_descriptions = load_pickle("smpdb_pathway_descriptions.pkl") or {}
        self.cached_uberon_mappings = load_pickle("uberon_tissue_name2id.pkl") or {}

    def _flatten_protein_functions(self, protein_function_list: List[Dict]) -> Dict:
        """Flatten list of protein function dictionaries into single dict."""
        flattened = {}
        for func_dict in protein_function_list:
            flattened.update(func_dict)
        return flattened

    def get_microbes(self, metabolite) -> List[str]:
        """Extract microbe names from metabolite ontology."""
        ontology = metabolite.find("hmdb:ontology", self.namespace)
        if ontology is None:
            return []

        for root in ontology.findall("hmdb:root", self.namespace):
            if self.get_text(root, "term") != "Disposition":
                continue

            for descendant in root.findall(".//hmdb:descendant", self.namespace):
                if self.get_text(descendant, "term") == "Microbe":
                    microbe_terms = descendant.findall(".//hmdb:term", self.namespace)
                    return sorted(
                        [
                            t.text.lower().strip()
                            for t in microbe_terms[1:]  # Skip first term == "Microbe"
                            if t.text and t.text.strip().lower() != "microbe"
                        ]
                    )
        return []

    def get_diseases(self, metabolite) -> List[str]:
        """Extract disease names from metabolite."""
        disease_names = set()
        diseases_elem = metabolite.find("hmdb:diseases", self.namespace)
        if diseases_elem is not None:
            for disease_elem in diseases_elem.findall("hmdb:disease", self.namespace):
                name_elem = disease_elem.find("hmdb:name", self.namespace)
                if name_elem is not None and name_elem.text:
                    disease_names.add(name_elem.text.strip().lower())
        return sorted(disease_names)

    def get_anatomical_entities(self, metabolite) -> Dict[str, List[str]]:
        """Extract anatomical entity information."""
        output = {
            "cellular_component": [],
            "biosample": [],
            "anatomical_entity": [],
        }

        bio_prop = metabolite.find("hmdb:biological_properties", self.namespace)
        if bio_prop is not None:
            # cellular locations
            cell_elem = bio_prop.find("hmdb:cellular_locations", self.namespace)
            if cell_elem is not None:
                for cell in cell_elem.findall("hmdb:cellular", self.namespace):
                    if cell.text and cell.text.strip():
                        output["cellular_component"].append(cell.text.strip().lower())

            # biospecimen locations
            specimen_elem = bio_prop.find("hmdb:biospecimen_locations", self.namespace)
            if specimen_elem is not None:
                for specimen in specimen_elem.findall("hmdb:biospecimen", self.namespace):
                    if specimen.text and specimen.text.strip():
                        specimen_name = specimen.text.strip().lower()
                        if specimen_name in self.cached_uberon_mappings:
                            output["biosample"].append(self.cached_uberon_mappings[specimen_name])
                        else:
                            output["biosample"].append(
                                {
                                    "original_name": specimen_name,
                                    "category": "biolink:AnatomicalEntity",
                                }
                            )

            # tissue locations
            tissue_elem = bio_prop.find("hmdb:tissue_locations", self.namespace)
            if tissue_elem is not None:
                for tissue in tissue_elem.findall("hmdb:tissue", self.namespace):
                    if tissue.text and tissue.text.strip():
                        tissue_name = tissue.text.strip().lower()
                        if tissue_name in self.cached_uberon_mappings:
                            output["anatomical_entity"].append(
                                self.cached_uberon_mappings[tissue_name]
                            )
                        else:
                            output["anatomical_entity"].append(
                                {
                                    "original_name": tissue_name,
                                    "category": "biolink:AnatomicalEntity",
                                }
                            )

        return output

    def get_publications(self, disease_elem) -> Dict:
        """Extract medical references from disease element."""
        pmids = []
        references_elem = disease_elem.find("hmdb:references", self.namespace)
        if references_elem is not None:
            for ref in references_elem.findall("hmdb:reference", self.namespace):
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

    def build_metabolite_node(
            self,
            metabolite,
            primary_id: str,
            xrefs: Dict,
            anatomical_entities: Dict = None,
    ) -> Dict:
        """Build standardized metabolite node."""
        name = self.get_text(metabolite, "name")
        state = self.get_text(metabolite, "state")

        synonyms_elem = metabolite.find("hmdb:synonyms", self.namespace)
        synonyms = self.get_list(synonyms_elem, "synonym") if synonyms_elem is not None else []

        node = {
            "id": primary_id,
            "name": name.lower() if name else None,
            "synonyms": synonyms,
            "description": self.get_text(metabolite, "description"),
            "chemical_formula": self.get_text(metabolite, "chemical_formula"),
            "molecular_weight": MolecularWeightExtractor.get_molecular_weights(metabolite, self),
            "state": state.lower() if state else None,
            "water_solubility": PropertyExtractor.get_experimental_property(
                metabolite, self, "water_solubility"
            ),
            "logp": PropertyExtractor.get_experimental_property(metabolite, self, "logp"),
            "melting_point": PropertyExtractor.get_experimental_property(
                metabolite, self, "melting_point"
            ),
            "category": "biolink:SmallMolecule",
            "xrefs": xrefs,
        }

        if anatomical_entities:
            node.update(anatomical_entities)

        return self.remove_empty_none_values(node)

    def parse_microbe_metabolite(self) -> Iterator[Dict]:
        """Parse microbe-metabolite associations."""

        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}metabolite':
                metabolite = elem

                microbes = self.get_microbes(metabolite)
                if not microbes:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                description = self.get_text(metabolite, "description")
                references = self.reference_extractor.extract_references(description, microbes)

                association_node = {
                    "category": "biolink:OrganismTaxonToChemicalEntityAssociation",  # Need to add this to Biolink Model
                    "predicate": "biolink:has_metabolic_interaction_with",  # Need to add this to Biolink Model
                    "primary_knowledge_source": "infores:hmdb_v5",
                    "has_evidence": "ECO:0000305" if references else "ECO:0000000",
                    "agent_type": "biolink:manual_agent",
                    "publication": references if references else None,
                }
                association_node = self.remove_empty_none_values(association_node)

                # metabolite node
                primary_id, xrefs = IDHierarchy.get_metabolite_primary_id(metabolite, self)
                object_node = self.build_metabolite_node(
                    metabolite,
                    primary_id,
                    xrefs,
                    anatomical_entities=self.get_anatomical_entities(metabolite),
                )

                # generate associations for each microbe
                for microbe in microbes:
                    if microbe in self.cached_taxon_info:
                        subject_node = self.cached_taxon_info[microbe].copy()
                        subject_node["original_name"] = microbe.lower().strip()
                        subject_node["category"] = "biolink:OrganismTaxon"
                        subject_node["organism_type"] = OrganismClassifier.get_organism_type(
                            subject_node
                        )
                        subject_node = self.remove_empty_none_values(subject_node)

                        _id = (
                            f"{subject_node['id'].split(':')[1]}"
                            f"_{association_node['predicate'].split(':')[1]}"
                            f"_{object_node['id'].split(':')[1]}"
                            if object_node.get("id", f"uuid:{str(uuid.uuid4())}")
                               and subject_node.get("id", f"uuid:{str(uuid.uuid4())}")
                            else f"uuid:{str(uuid.uuid4())}"
                        )

                        yield {
                            "_id": _id,
                            "association": association_node,
                            "object": object_node,
                            "subject": subject_node,
                        }

                # clear memory after processing all metabolites
                elem.clear()
                for ancestor in elem.xpath("ancestor-or-self::*"):
                    ancestor.clear()

    def parse_metabolite_disease(self) -> Iterator[Dict]:
        """Parse metabolite-disease associations."""
        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}metabolite':
                metabolite = elem

                # metabolite node
                primary_id, xrefs = IDHierarchy.get_metabolite_primary_id(metabolite, self)
                subject_node = self.build_metabolite_node(
                    metabolite,
                    primary_id,
                    xrefs,
                    anatomical_entities=self.get_anatomical_entities(metabolite),
                )

                diseases = self.get_diseases(metabolite)
                if not diseases:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                diseases_elem = metabolite.find("hmdb:diseases", self.namespace)
                if diseases_elem is None:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                for disease_elem in diseases_elem.findall("hmdb:disease", self.namespace):
                    references = self.get_publications(disease_elem)

                    # association node
                    association_node = {
                        "category": "biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation",
                        "predicate": "biolink:associated_with",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "has_evidence": "ECO:0000305" if references else "ECO:0000000",
                        "agent_type": "biolink:manual_agent",
                        "publication": references if references else None,
                    }
                    association_node = self.remove_empty_none_values(association_node)

                    # disease node
                    for disease in diseases:
                        disease_key = disease.strip().lower()
                        if disease_key in self.cached_disease_info:
                            object_node = self.cached_disease_info[disease_key].copy()
                            object_node = self.remove_empty_none_values(object_node)

                            _id = (
                                f"{subject_node['id'].split(':')[1]}"
                                f"_{association_node['predicate'].split(':')[1]}"
                                f"_{object_node['id'].split(':')[1]}"
                                if object_node.get("id", f"uuid:{str(uuid.uuid4())}")
                                   and subject_node.get("id", f"uuid:{str(uuid.uuid4())}")
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

    def parse_metabolite_protein(self) -> Iterator[Dict]:
        """Parse metabolite-protein associations."""

        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}metabolite':
                metabolite = elem
                # metabolite node
                primary_id, xrefs = IDHierarchy.get_metabolite_primary_id(metabolite, self)
                subject_node = self.build_metabolite_node(
                    metabolite,
                    primary_id,
                    xrefs,
                    anatomical_entities=self.get_anatomical_entities(metabolite),
                )

                # process protein associations
                prot_elem = metabolite.find("hmdb:protein_associations", self.namespace)
                if prot_elem is None:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                for protein in prot_elem.findall("hmdb:protein", self.namespace):
                    # protein node
                    object_node = {
                        "id": None,
                        "name": self.get_text(protein, "gene_name"),
                        "full_name": self.get_text(protein, "name"),
                        "description": None,
                        "category": "biolink:Protein",
                        "protein_type": None,
                        "xrefs": {},
                    }

                    # add protein type
                    prot_type = self.get_text(protein, "protein_type")
                    if prot_type:
                        object_node["protein_type"] = prot_type.lower()

                    # add HMDBP accession
                    prot_accession = self.get_text(protein, "protein_accession")
                    if prot_accession:
                        object_node["xrefs"]["hmdbp"] = f"HMDBP:{prot_accession}"

                    # UniProt ID and description
                    uniprot_id = self.get_text(protein, "uniprot_id")
                    if uniprot_id:
                        object_node["id"] = f"UniProtKB:{uniprot_id}"
                        object_node["xrefs"]["uniprotkb"] = f"UniProtKB:{uniprot_id}"

                        # add protein function if available
                        if uniprot_id in self.cached_protein_functions:
                            object_node["description"] = self.cached_protein_functions[uniprot_id][
                                "description"
                            ]
                    elif prot_accession:
                        object_node["id"] = f"HMDBP:{prot_accession}"

                    object_node = self.remove_empty_none_values(object_node)

                    # association node
                    association_node = {
                        "category": "biolink:ChemicalGeneInteractionAssociation",
                        "predicate": "biolink:interacts_with",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "has_evidence": "ECO:0000000",  # unknown evidence
                        "agent_type": "biolink:manual_agent",
                    }

                    _id = (
                        f"{subject_node['id'].split(':')[1]}"
                        f"_{association_node['predicate'].split(':')[1]}"
                        f"_{object_node['id'].split(':')[1]}"
                        if object_node.get("id", f"uuid:{str(uuid.uuid4())}")
                           and subject_node.get("id", f"uuid:{str(uuid.uuid4())}")
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

    def parse_metabolite_pathway(self) -> Iterator[Dict]:
        """Parse metabolite-pathway associations."""
        for event, elem in ET.iterparse(self.input_xml, events=("start", "end")):
            if event == "end" and elem.tag == f'{{{self.namespace["hmdb"]}}}metabolite':
                metabolite = elem

                # metabolite node
                primary_id, xrefs = IDHierarchy.get_metabolite_primary_id(metabolite, self)
                subject_node = self.build_metabolite_node(
                    metabolite,
                    primary_id,
                    xrefs,
                    anatomical_entities=self.get_anatomical_entities(metabolite),
                )

                # pathway associations
                bio_prop = metabolite.find("hmdb:biological_properties", self.namespace)
                if bio_prop is None:
                    elem.clear()
                    for ancestor in elem.xpath("ancestor-or-self::*"):
                        ancestor.clear()
                    continue

                pathways_elem = bio_prop.find("hmdb:pathways", self.namespace)
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

                    # get SMPDB info from cache
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
                        "category": "biolink:ChemicalToPathwayAssociation",
                        "predicate": "biolink:participates_in",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "has_evidence": "ECO:0000000",  # unknown evidence
                        "agent_type": "biolink:manual_agent",
                    }

                    # generate association ID
                    _id = (
                        f"{subject_node['id'].split(':')[1]}"
                        f"_{association_node['predicate'].split(':')[1]}"
                        f"_{object_node['id'].split(':')[1]}"
                        if object_node.get("id", f"uuid:{str(uuid.uuid4())}")
                           and subject_node.get("id", f"uuid:{str(uuid.uuid4())}")
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
