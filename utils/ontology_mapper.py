"""Centralized ontology mapping orchestration with fallback strategies."""

import asyncio
from typing import Dict, List

from manual_annotations.disease_name2id import format_manual_disease_mappings
from manual_annotations.taxon_name2taxid import (
    apply_text2term_corrections,
    format_manual_mappings,
)
from .cache_helper import load_pickle, save_pickle
from .ontology_services import (
    BiothingsServices,
    DiseaseServices,
    NCITServices,
    TaxonServices,
    UMLSClient,
)


class TaxonMapper:
    """Orchestrates taxon name to taxonomy ID mapping with multiple strategies."""

    def __init__(self, email: str):
        self.email = email

    def map_all_taxon_names(self, taxon_names: List[str], use_cache: bool = True) -> Dict:
        """Map taxon names using multiple strategies with fallbacks."""
        if use_cache:
            cached = load_pickle("all_taxon_name2taxid.pkl")
            if cached:
                return cached

        # step 1: ETE3 (fast, exact matches)
        print("Step 1: ETE3 taxon mapping...")
        ete3_mapped = TaxonServices.ete3_taxon_name2taxid(taxon_names)
        save_pickle(ete3_mapped, "ete3_name2taxid.pkl")

        # step 2: Entrez (for ETE3 misses)
        no_hits = [name for name in taxon_names if name not in ete3_mapped]
        print(f"Step 2: Entrez mapping for {len(no_hits)} unmapped names...")
        entrez_mapped = TaxonServices.entrez_taxon_name2taxid(no_hits, self.email)
        save_pickle(entrez_mapped, "entrez_name2taxid.pkl")

        # step 3: Text2Term (for remaining misses)
        no_hits2 = [
            name for name in taxon_names if name not in ete3_mapped and name not in entrez_mapped
        ]
        print(f"Step 3: Text2term mapping for {len(no_hits2)} unmapped names...")
        text2term_mapped = TaxonServices.text2term_taxon_name2taxid(no_hits2)
        text2term_mapped = apply_text2term_corrections(text2term_mapped)
        save_pickle(text2term_mapped, "text2term_name2taxid.pkl")

        # step 4: manual mappings (for final misses)
        no_hits3 = [
            name
            for name in taxon_names
            if name not in ete3_mapped
            and name not in entrez_mapped
            and name not in text2term_mapped
        ]
        print(f"Step 4: Manual mapping for {len(no_hits3)} unmapped names...")
        manual_mapped = format_manual_mappings(no_hits3)
        save_pickle(manual_mapped, "manual_name2taxid.pkl")

        # combine all mappings
        all_mapped = ete3_mapped | entrez_mapped | text2term_mapped | manual_mapped
        save_pickle(all_mapped, "all_taxon_name2taxid.pkl")

        return all_mapped

    def enrich_taxon_info(self, mapped_taxons: Dict) -> Dict:
        """Enrich taxon mappings with detailed information."""
        # step 1: get taxonomy IDs
        taxids = [taxon_data["taxid"] for taxon_data in mapped_taxons.values()]

        # step 2: get detailed taxon info from biothings
        print("Enriching with detailed taxon information...")
        taxon_info = BiothingsServices.get_taxon_info_from_bt(taxids)

        # step 3: get scientific names for NCIT description lookup
        sci_names = [info["name"] for info in taxon_info.values() if "name" in info]

        # step 4: get NCIT descriptions
        print("Fetching NCIT descriptions...")
        taxon_descr = NCITServices.get_ncit_taxon_description(sci_names)

        # step5: add descriptions to taxon info
        for _taxid_str, info in taxon_info.items():
            name = info.get("name")
            if name and name in taxon_descr:
                info.update(taxon_descr[name])

        # step 6: create final mapping: original_name -> enriched_info
        enriched_taxon_info = {}
        for original_name, mapping_data in mapped_taxons.items():
            taxid_str = str(mapping_data["taxid"])
            if taxid_str in taxon_info:
                enriched_info = taxon_info[taxid_str].copy()
                enriched_info["original_name"] = original_name
                enriched_info["mapping_tool"] = mapping_data["mapping_tool"]
                enriched_taxon_info[original_name] = enriched_info

        save_pickle(enriched_taxon_info, "original_taxon_name2taxid.pkl")
        return enriched_taxon_info


class DiseaseMapper:
    """Orchestrates disease name to ID mapping with multiple strategies."""

    def __init__(self, umls_api_key: str):
        self.umls_client = UMLSClient(umls_api_key)

    def map_all_disease_names(self, disease_dict: Dict[str, int], use_cache: bool = True) -> Dict:
        """Map disease names using multiple strategies with fallbacks."""
        if use_cache:
            cached = load_pickle("all_disease_name2id.pkl")
            if cached:
                return cached

        # step 1: HMDB OMIM mappings (direct from data)
        print("Step 1: Processing HMDB OMIM mappings...")
        hmdb_mapped = {
            name: {"id": f"OMIM:{omim}", "mapping_tool": "hmdb_v5"}
            for name, omim in disease_dict.items()
            if omim
        }

        # step 2: Text2Term for unmapped diseases
        unmapped_diseases = [name for name, omim in disease_dict.items() if not omim]
        print(f"Step 2: Text2term mapping for {len(unmapped_diseases)} diseases...")
        text2term_mapped = DiseaseServices.text2term_disease_name2id(unmapped_diseases)
        save_pickle(text2term_mapped, "text2term_disease_name2id.pkl")

        # step 3: UMLS for Text2Term misses
        umls_candidates = [name for name in unmapped_diseases if name not in text2term_mapped]
        print(f"Step 3: UMLS mapping for {len(umls_candidates)} diseases...")
        umls_mapped = asyncio.run(self.umls_client.query_cuis(umls_candidates))
        save_pickle(umls_mapped, "umls_disease_name2id.pkl")

        # step 4: manual mappings for final misses
        manual_candidates = [name for name, mapped in umls_mapped.items() if not mapped.get("id")]
        print(f"Step 4: Manual mapping for {len(manual_candidates)} diseases...")
        manual_mapped = format_manual_disease_mappings(manual_candidates)
        save_pickle(manual_mapped, "manual_disease_name2id.pkl")

        # step 5: combine all mappings, filter out empty IDs
        all_mapped = hmdb_mapped | text2term_mapped | umls_mapped | manual_mapped
        filtered_mapped = {
            name: mapped
            for name, mapped in all_mapped.items()
            if mapped.get("id") and mapped["id"] != ""
        }

        save_pickle(filtered_mapped, "all_disease_name2id.pkl")
        return filtered_mapped

    def enrich_disease_info(self, mapped_diseases: Dict) -> Dict:
        """Enrich disease mappings with detailed information."""
        # step 1: extract IDs for biothings lookup
        bt_ids = []
        for mapped in mapped_diseases.values():
            _id = mapped["id"]
            if "UMLS" not in _id and "OMIM" not in _id:
                bt_ids.append(_id)
            else:
                bt_ids.append(_id.split(":")[1])

        # step 2: get detailed disease info from biothings
        print("Enriching with detailed disease information...")
        bt_disease_info = BiothingsServices.bt_get_disease_info(bt_ids)

        # step 3: create enriched mapping
        enriched_disease_info = {}
        for original_name, mapping_data in mapped_diseases.items():
            _id = mapping_data["id"]
            if _id in bt_disease_info:
                info = bt_disease_info[_id].copy()
                info["original_name"] = original_name
                info["mapping_tool"] = mapping_data["mapping_tool"]
                enriched_disease_info[original_name] = info
            else:
                # create basic info for unmapped diseases
                enriched_disease_info[original_name] = {
                    "id": _id,
                    "name": original_name,
                    "type": "biolink:Disease",
                    "original_name": original_name,
                    "mapping_tool": mapping_data["mapping_tool"],
                }

        save_pickle(enriched_disease_info, "original_disease_name2id.pkl")
        return enriched_disease_info


class ProteinMapper:
    """Orchestrates protein and gene mapping with enrichment."""

    @staticmethod
    def enrich_protein_mappings(protein2uniprot: Dict[str, str]) -> Dict:
        """Enrich protein mappings with function descriptions."""
        from .ontology_services import ProteinServices

        # step 1: get UniProt IDs
        uniprot_ids = [uid for uid in protein2uniprot.values() if uid]

        # step 2/1: get protein functions
        print("Fetching protein functions...")
        protein_functions = asyncio.run(ProteinServices.get_batch_protein_functions(uniprot_ids))
        save_pickle(protein_functions, "uniprot_protein_functions.pkl")

        # step 2/2: get gene mappings
        print("Mapping UniProt to EntrezGene...")
        uniprot2entrez = ProteinServices.uniprot_id2entrezgene(uniprot_ids)
        save_pickle(uniprot2entrez, "bt_uniprot2entrezgene.pkl")

        return {"protein_functions": protein_functions, "uniprot2entrez": uniprot2entrez}

    @staticmethod
    def enrich_gene_info(uniprot2entrez: Dict) -> Dict:
        """Enrich gene mappings with summaries."""
        from .ontology_services import GeneServices

        # step 1: extract Entrez gene IDs
        entrez_genes = [
            mapped["gene_id"].split(":")[1]
            for mapped in uniprot2entrez.values()
            if "gene_id" in mapped
        ]

        # step 2: get gene summaries
        print("Fetching gene summaries...")
        gene_summaries = asyncio.run(GeneServices.get_batch_gene_summaries(entrez_genes))
        save_pickle(gene_summaries, "entrezgene_summaries.pkl")

        return gene_summaries
