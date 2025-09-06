"""
Cache management for HMDB data processing pipeline.
"""
import asyncio
import os
import sys
from typing import Dict, List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))

from manual_annotations.anatomical_name2uberonid import apply_manual_anatomy_mappings

from .cache_helper import cache_exists, load_pickle, save_pickle
from .ontology_mapper import DiseaseMapper, ProteinMapper, TaxonMapper
from .ontology_services import (
    GeneServices,
    GOServices,
    PathwayServices,
    ProteinServices,
    UberonService,
)
from .reader import (
    get_all_anatomical_terms_from_hmdb,
    get_all_diseases,
    get_all_go_terms_from_hmdbp,
    get_all_microbe_names,
    get_all_uniprot_ids_from_hmdb,
    get_all_uniprot_ids_from_hmdbp,
)


class CacheManager:
    """Manages all caching operations for HMDB data processing."""

    def __init__(self, email: str, umls_api_key: str):
        self.email = email
        self.umls_api_key = umls_api_key
        self.taxon_mapper = TaxonMapper(email)
        self.disease_mapper = DiseaseMapper(umls_api_key)

    def cache_metabolite_data(self, metabolite_xml: str, skip_existing: bool = True) -> None:
        """Cache all metabolite-related mappings and enrichments with selective skipping."""
        print("\n>>> Caching Metabolite Data...")
        print("=" * 50)

        # step 1: cache taxon mappings
        if not skip_existing or not cache_exists("original_taxon_name2taxid.pkl"):
            print(">>> Processing microbe taxon names...")
            microbe_names = list(set(get_all_microbe_names(metabolite_xml)))
            save_pickle(microbe_names, "hmdb_v5_microbe_names.pkl")

            # map taxon names using ete3 -> entrez -> text2term -> manual
            all_mapped_taxons = self.taxon_mapper.map_all_taxon_names(microbe_names)

            # enrich with detailed taxon info
            enriched_taxon_info = self.taxon_mapper.enrich_taxon_info(all_mapped_taxons)
            print(f"[DONE] Mapped {len(enriched_taxon_info)} taxon names")
        else:
            print("[DONE] Taxon mappings already cached, skipping...")

        # step 2: cache disease mappings
        if not skip_existing or not cache_exists("original_disease_name2id.pkl"):
            print(">>>️ Processing disease names...")
            diseases = get_all_diseases(metabolite_xml)
            save_pickle(diseases, "hmdb_v5_diseases.pkl")

            # map disease names using text2term -> umls -> manual
            all_mapped_diseases = self.disease_mapper.map_all_disease_names(diseases)

            # enrich with detailed disease info
            enriched_disease_info = self.disease_mapper.enrich_disease_info(all_mapped_diseases)
            print(f"[DONE] Mapped {len(enriched_disease_info)} disease names")
        else:
            print("[DONE] Disease mappings already cached, skipping...")

        # step 3: cache protein mappings
        if not skip_existing or not cache_exists("uniprot_protein_functions.pkl"):
            print(">>>️ Processing protein mappings...")
            protein2uniprot = get_all_uniprot_ids_from_hmdb(metabolite_xml)
            save_pickle(protein2uniprot, "all_protein_name2uniprot.pkl")

            # enrich protein mappings
            protein_enrichment = ProteinMapper.enrich_protein_mappings(protein2uniprot)

            # enrich gene information
            ProteinMapper.enrich_gene_info(protein_enrichment["uniprot2entrez"])
            print(f"[DONE] Processed {len(protein2uniprot)} protein mappings")
        else:
            print("[DONE] Protein mappings already cached, skipping...")

        # step 4: cache pathway descriptions
        if not skip_existing or not cache_exists("smpdb_pathway_descriptions.pkl"):
            print("\n>>> Processing pathway descriptions...")
            pathway_descriptions = PathwayServices.get_smpdb_pathway_description()
            save_pickle(pathway_descriptions, "smpdb_pathway_descriptions.pkl")
            print(f"[DONE] Cached {len(pathway_descriptions)} pathway descriptions")
        else:
            print("[DONE] Pathway descriptions already cached, skipping...")

        # step 5: cache Uberon anatomy terms (from tissue name to Uberon ID)
        if not skip_existing or not cache_exists("uberon_tissue_name2id.pkl"):
            print(">>>️ Querying UBERON anatomy terms...")
            anatomical_terms = get_all_anatomical_terms_from_hmdb(metabolite_xml)
            if not anatomical_terms:
                print("!! No anatomical terms found in metabolite data")
                return
            else:
                print(f"-> Found {len(anatomical_terms)} unique anatomical terms")
            uberon_terms = UberonService.query_terms(anatomical_terms)
            corrected_uberon_terms = apply_manual_anatomy_mappings(uberon_terms)
            save_pickle(corrected_uberon_terms, "uberon_tissue_name2id.pkl")
            print(f"-> Cached {len(corrected_uberon_terms)} UBERON anatomy terms")
        else:
            print("[DONE] UBERON anatomy terms already cached, skipping...")

        print("\n[DONE] Metabolite Data Caching Complete!")

    def cache_protein_data(self, protein_xml: str, skip_existing: bool = True) -> None:
        """Cache all protein-related mappings and enrichments with selective skipping."""
        print("\n>>> Caching Protein Data...")
        print("=" * 50)

        # step 1: cache HMDBP protein functions
        if not skip_existing or not cache_exists("hmdbp_uniprot_protein_functions.pkl"):
            print(">>>️ Processing HMDBP UniProt IDs...")
            hmdbp_uniprot_ids = get_all_uniprot_ids_from_hmdbp(protein_xml)

            print(">>> Querying HMDBP protein functions...")
            hmdbp_protein_functions = asyncio.run(
                ProteinServices.get_batch_protein_functions(hmdbp_uniprot_ids)
            )
            save_pickle(hmdbp_protein_functions, "hmdbp_uniprot_protein_functions.pkl")
        else:
            print("[DONE] HMDBP protein functions already cached, skipping...")

        # step 2: cache gene mappings and summaries
        if not skip_existing or not cache_exists("hmdbp_uniprot2entrezgene.pkl"):
            print(">>> Mapping HMDBP UniProt to EntrezGene...")
            hmdbp_uniprot_ids = get_all_uniprot_ids_from_hmdbp(protein_xml)
            hmdbp_uniprot2entrez = ProteinMapper.uniprot_id2entrezgene(hmdbp_uniprot_ids)
            save_pickle(hmdbp_uniprot2entrez, "hmdbp_uniprot2entrezgene.pkl")
        else:
            print("[DONE] HMDBP UniProt to EntrezGene mappings already cached, skipping...")
            hmdbp_uniprot2entrez = load_pickle("hmdbp_uniprot2entrezgene.pkl")

        if not skip_existing or not cache_exists("hmdbp_entrezgene_summaries.pkl"):
            # get gene summaries
            entrez_genes = [
                gene_info["gene_id"].split(":")[1]
                for gene_info in hmdbp_uniprot2entrez.values()
                if "gene_id" in gene_info
            ]

            print(">>> Querying HMDBP gene summaries...")
            hmdbp_gene_summaries = asyncio.run(GeneServices.get_batch_gene_summaries(entrez_genes))
            save_pickle(hmdbp_gene_summaries, "hmdbp_entrezgene_summaries.pkl")
        else:
            print("[DONE] HMDBP gene summaries already cached, skipping...")

        # step 3: cache GO terms and definitions
        if not skip_existing or not cache_exists("hmdbp_go_definitions.pkl"):
            print(">>> Processing GO terms...")
            go_terms = get_all_go_terms_from_hmdbp(protein_xml)

            print(">>> Querying GO definitions...")
            go_definitions = asyncio.run(GOServices.get_go_definitions(go_terms))
            save_pickle(go_definitions, "hmdbp_go_definitions.pkl")
            print(f"[DONE] Cached {len(go_definitions)} GO definitions")
        else:
            print("[DONE] GO definitions already cached, skipping...")

        print("\n[DONE] Protein Data Caching Complete!")

    def get_cache_status(self) -> Dict[str, bool]:
        """Check which cache files exist."""
        cache_files = [
            "hmdb_v5_microbe_names.pkl",
            "ete3_name2taxid.pkl",
            "entrez_name2taxid.pkl",
            "text2term_name2taxid.pkl",
            "manual_name2taxid.pkl",
            "all_taxon_name2taxid.pkl",
            "original_taxon_name2taxid.pkl",
            "hmdb_v5_diseases.pkl",
            "text2term_disease_name2id.pkl",
            "umls_disease_name2id.pkl",
            "manual_disease_name2id.pkl",
            "all_disease_name2id.pkl",
            "original_disease_name2id.pkl",
            "all_protein_name2uniprot.pkl",
            "uniprot_protein_functions.pkl",
            "bt_uniprot2entrezgene.pkl",
            "entrezgene_summaries.pkl",
            "smpdb_pathway_descriptions.pkl",
            "hmdbp_uniprot_protein_functions.pkl",
            "hmdbp_uniprot2entrezgene.pkl",
            "hmdbp_entrezgene_summaries.pkl",
            "hmdbp_go_definitions.pkl",
            "uberon_tissue_name2id.pkl",
        ]

        return {filename: cache_exists(filename) for filename in cache_files}

    def is_cache_complete(self) -> bool:
        """Check if all required cache files exist."""
        status = self.get_cache_status()
        return all(status.values())

    def get_missing_cache_files(self) -> List[str]:
        """Get list of missing cache files."""
        status = self.get_cache_status()
        return [filename for filename, exists in status.items() if not exists]

    def load_cached_data(self) -> Dict:
        """Load all cached data for parsing."""
        return {
            "taxon_info": load_pickle("original_taxon_name2taxid.pkl"),
            "disease_info": load_pickle("original_disease_name2id.pkl"),
            "protein_functions": load_pickle("uniprot_protein_functions.pkl"),
            "pathway_descriptions": load_pickle("smpdb_pathway_descriptions.pkl"),
            "hmdbp_protein_functions": load_pickle("hmdbp_uniprot_protein_functions.pkl"),
            "hmdbp_uniprot2entrez": load_pickle("hmdbp_uniprot2entrezgene.pkl"),
            "hmdbp_gene_summaries": load_pickle("hmdbp_entrezgene_summaries.pkl"),
            "go_definitions": load_pickle("hmdbp_go_definitions.pkl"),
            "uberon_terms": load_pickle("uberon_tissue_name2id.pkl"),
        }

    def cache_exists_for_file(self, xml_file: str) -> bool:
        """Check if cache exists for a specific XML file type."""
        if "metabolite" in xml_file.lower():
            required_files = [
                "original_taxon_name2taxid.pkl",
                "original_disease_name2id.pkl",
                "uniprot_protein_functions.pkl",
                "smpdb_pathway_descriptions.pkl",
                "uberon_tissue_name2id.pkl",
            ]
        elif "protein" in xml_file.lower():
            required_files = ["hmdbp_uniprot_protein_functions.pkl", "hmdbp_go_definitions.pkl"]
        else:
            return False

        return all(cache_exists(f) for f in required_files)
