"""Helper utilities for XML parsing and data transformation."""

import re
from typing import Any, Dict, Union


class XMLParseHelper:
    """Base helper class for XML parsing operations."""

    namespace = {"hmdb": "http://www.hmdb.ca"}

    def __init__(self, input_xml: str):
        self.input_xml = input_xml

    def get_text(self, elem, tag: str) -> Union[str, None]:
        """Extract text from XML element."""
        child = elem.find(f"hmdb:{tag}", self.namespace)
        return child.text.strip() if child is not None and child.text else None

    def get_list(self, elem, tag: str) -> list[str]:
        """Extract list of text values from XML elements."""
        return [e.text.lower() for e in elem.findall(f"hmdb:{tag}", self.namespace) if e.text]

    def remove_empty_none_values(self, obj: Any) -> Any:
        """Recursively remove empty/None values from dictionaries and lists."""
        if isinstance(obj, dict):
            cleaned = {}
            for k, v in obj.items():
                v_clean = self.remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned[k] = v_clean
            return cleaned

        if isinstance(obj, list):
            cleaned_list = []
            for v in obj:
                v_clean = self.remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned_list.append(v_clean)
            return cleaned_list
        return obj


class OrganismClassifier:
    """Utilities for organism classification."""

    @staticmethod
    def get_organism_type(node: Dict) -> str:
        """Classify the organism type based on lineage taxonomy IDs."""
        taxon_map = {
            2: "biolink:Bacterium",
            2157: "Archaeon",
            2759: "Eukaryote",
            10239: "biolink:Virus",
        }

        for taxid, biolink_type in taxon_map.items():
            if taxid in node.get("lineage", []):
                return biolink_type

        return "Other"


class ReferenceExtractor:
    """Extract and format reference information from the text."""

    def __init__(self):
        self.parenthetical_pattern = re.compile(r"([^.?!]*?)\s*\(([^)]*?)\)")
        self.reference_indicators = ["PMID", "DOI", "Wikipedia", "http", "www", ":"]

    def extract_references(self, description: str, microbes: list = None) -> Dict:
        """Extract structured reference information from description text."""
        if not description:
            return {}

        microbes = microbes or []
        matches = self.parenthetical_pattern.findall(description)

        for sentence, paren in matches:
            if not any(ind.lower() in paren.lower() for ind in self.reference_indicators):
                continue

            # check if any microbes are mentioned in this sentence
            matched_microbes = [m for m in microbes if m.lower() in sentence.lower()]
            if not matched_microbes:
                continue

            # clean up PMID references
            paren = re.sub(r"\bPMID[\s:]+(\d+)", r"PMID:\1", paren, flags=re.IGNORECASE)
            refs = [r.strip() for r in re.split(r"[;|,]", paren) if ":" in r and "CAS" not in r]

            if refs:
                return self._format_references(refs)

        return {}

    def _format_references(self, refs: list) -> Dict:
        """Format the list of references into a structured dictionary."""
        ref_dict = {}

        for ref in refs:
            try:
                prefix, value = ref.split(":", 1)
            except ValueError:
                continue

            key = self._classify_reference_type(prefix)
            value = self._clean_reference_value(key, value)

            ref_dict.setdefault(key, []).append(value)

        # convert single-item lists to scalars
        for k in list(ref_dict):
            if isinstance(ref_dict[k], list) and len(ref_dict[k]) == 1:
                ref_dict[k] = ref_dict[k][0]

        ref_dict["category"] = "biolink:Publication"

        return ref_dict

    def _classify_reference_type(self, prefix: str) -> str:
        """Classify reference prefix into standard types."""
        prefix_lower = prefix.lower()

        if prefix_lower == "pmid":
            return "pmid"
        elif prefix_lower == "doi":
            return "doi"
        elif "wikipedia" in prefix_lower:
            return "wikidata"
        elif "http" in prefix_lower or "www" in prefix_lower or ".com" in prefix_lower:
            return "url"
        else:
            return "article"

    def _clean_reference_value(self, ref_type: str, value: str) -> Union[str, int]:
        """Clean and format reference values."""
        value = value.strip()

        if ref_type == "pmid" and value.isdigit():
            return f"PMID:{int(value)}"
        elif ref_type == "doi" and ":" in value:
            return value.strip()

        return value

    def _generate_primary_id(self, ref_dict: Dict) -> str:
        """Generate primary ID based on reference priority."""
        if "pmid" in ref_dict:
            pmid = ref_dict["pmid"]
            if isinstance(pmid, list):
                pmid = pmid[0]
            return f"PMID:{pmid}"
        elif "doi" in ref_dict:
            doi = ref_dict["doi"]
            if isinstance(doi, list):
                doi = doi[0]
            return f"doi:{doi}"
        elif "url" in ref_dict:
            url = ref_dict["url"]
            if isinstance(url, list):
                url = url[0]
            return f"JournalArticle:{url}"
        elif "wikidata" in ref_dict:
            return "Wikipedia"
        elif "article" in ref_dict:
            article = ref_dict["article"]
            if isinstance(article, list):
                article = article[0]
            return str(article)

        return "Unknown"


class IDHierarchy:
    """Handle ID hierarchy and primary ID selection."""

    @staticmethod
    def get_metabolite_primary_id(metabolite, xml_helper: XMLParseHelper) -> tuple[str, Dict]:
        """Extract primary ID and xrefs for metabolite with hierarchy."""
        id_hierarchy = [
            ("pubchem_compound_id", "PUBCHEM.COMPOUND"),
            ("inchikey", "INCHIKEY"),
            ("smiles", None),
            ("drugbank_id", "DRUGBANK"),
            ("chebi_id", "CHEBI"),
            ("chembl_id", "CHEMBL.COMPOUND"),
            ("accession", "HMDB"),
            ("cas_registry_number", "CAS"),
            ("kegg_id", None),
            ("metlin_id", "METLIN"),
            ("chemspider_id", "chemspider"),
            ("foodb_id", "foodb.compound"),
            ("bigg_id", "BIGG.METABOLITE"),
            ("pdb_id", "PDB"),
            ("vmh_id", "VMH"),
        ]

        def classify_kegg(val: str) -> str:
            if val.startswith("C"):
                return "KEGG.COMPOUND"
            elif val.startswith("G"):
                return "KEGG.GLYCAN"
            elif val.startswith("D"):
                return "KEGG.DRUG"
            return "KEGG"

        xrefs, primary_id = {}, None

        for tag, prefix in id_hierarchy:
            val = xml_helper.get_text(metabolite, tag)
            if not val:
                continue

            if tag == "kegg_id":
                prefix = classify_kegg(val)
            elif tag == "smiles":
                prefix = ""

            curie = f"{prefix}:{val}"

            # determine xref key
            if prefix == "PUBCHEM.COMPOUND":
                key = "pubchem_cid"
            elif prefix == "HMDB":
                key = "hmdb"
            elif prefix == "METLIN":
                key = "metlin"
            else:
                key = prefix.lower().split(".")[0] if prefix else "smiles"

            if primary_id is None:
                primary_id = curie

            xrefs.setdefault(key, curie)

        return primary_id, xrefs

    @staticmethod
    def get_protein_primary_id(protein, xml_helper: XMLParseHelper) -> tuple[str, Dict]:
        """Extract primary ID and xrefs for protein with hierarchy."""
        id_hierarchy = [
            ("uniprot_id", "UniProtKB"),
            ("hgnc_id", "HGNC"),
            ("genbank_protein_id", "GBP"),
            ("genbank_gene_id", "GBG"),
            ("genecard_id", "GENECARD"),
            ("geneatlas_id", "GENEATLAS"),
            ("accession", "HMDBP"),
        ]

        xrefs, primary_id = {}, None

        for tag, prefix in id_hierarchy:
            val = xml_helper.get_text(protein, tag)
            if not val:
                continue

            curie = f"{prefix}:{val}"

            # determine xref key
            if prefix == "GBP":
                key = "genbank_prot"
            elif prefix == "HMDBP":
                key = "hmdbp"
            elif prefix == "GBG":
                key = "genbank_gene"
            else:
                key = prefix.lower()

            if primary_id is None:
                primary_id = curie

            xrefs.setdefault(key, curie)

        return primary_id, xrefs


class MolecularWeightExtractor:
    """Extract and format molecular weight information."""

    @staticmethod
    def get_molecular_weights(metabolite, xml_helper: XMLParseHelper) -> Dict[str, float]:
        """Extract average and monoisotopic molecular weights."""
        avg = xml_helper.get_text(metabolite, "average_molecular_weight")
        mono = xml_helper.get_text(metabolite, "monisotopic_molecular_weight")

        weights = {}
        if avg:
            try:
                weights["average_molecular_weight"] = float(avg)
            except ValueError:
                pass

        if mono:
            try:
                weights["monoisotopic_molecular_weight"] = float(mono)
            except ValueError:
                pass

        return weights


class PropertyExtractor:
    """Extract experimental and calculated properties."""

    @staticmethod
    def get_experimental_property(
            metabolite, xml_helper: XMLParseHelper, prop_name: str
    ) -> Union[str, None]:
        """Extract specific experimental property by name."""
        props = metabolite.find("hmdb:experimental_properties", xml_helper.namespace)
        if props is not None:
            for prop in props.findall("hmdb:property", xml_helper.namespace):
                kind = xml_helper.get_text(prop, "kind")
                if kind and kind.lower() == prop_name.lower():
                    return xml_helper.get_text(prop, "value")
        return None
