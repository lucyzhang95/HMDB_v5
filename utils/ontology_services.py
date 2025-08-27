"""Ontology and external API services with async support and progress tracking."""
import asyncio
import os
import ssl
import time
from typing import Dict, List, Optional

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
import text2term
from Bio import Entrez
from ete3 import NCBITaxa
from tqdm.auto import tqdm


class TaxonServices:
    """Services for taxonomy name to ID mapping."""

    @staticmethod
    def ete3_taxon_name2taxid(taxon_names: list) -> dict:
        """Map taxonomy names to NCBI taxonomy IDs using ete3."""
        ete3_mapped = {}
        taxon_names = set(taxon_names)

        if not NCBITaxa:
            ssl._create_default_https_context = ssl._create_unverified_context
            ncbi = NCBITaxa()
            ncbi.update_taxonomy_database()

        ncbi = NCBITaxa()

        with tqdm(desc="ETE3 taxon mapping", unit="names") as pbar:
            ete3_name2taxid = ncbi.get_name_translator(taxon_names)
            for name, taxid in ete3_name2taxid.items():
                if taxid:
                    ete3_mapped[name] = {"taxid": int(taxid[0]), "mapping_tool": "ete3"}
                pbar.update(1)

        return ete3_mapped

    @staticmethod
    def entrez_taxon_name2taxid(
        taxon_names: list[str],
        email: str,
        sleep: float = 0.34,
        retries: int = 3,
        backoff_factor: int = 2,
    ) -> dict:
        """Map taxonomy names to NCBI taxonomy IDs using Entrez API."""
        Entrez.email = email
        entrez_mapped = {}

        with tqdm(desc="Entrez taxon mapping", total=len(set(taxon_names)), unit="names") as pbar:
            for name in set(taxon_names):
                delay = 5
                for attempt in range(retries):
                    try:
                        handle = Entrez.esearch(db="taxonomy", term=name, retmode="xml", retmax=1)
                        record = Entrez.read(handle)
                        handle.close()
                        if record["IdList"]:
                            taxid = int(record["IdList"][0])
                            entrez_mapped[name] = {"taxid": taxid, "mapping_tool": "entrez"}
                        break
                    except Exception as e:
                        tqdm.write(
                            f"Entrez query failed for '{name}' on attempt {attempt + 1}/{retries}: {e}"
                        )
                        if attempt < retries - 1:
                            time.sleep(delay)
                            delay *= backoff_factor
                        else:
                            tqdm.write(f"All retry attempts failed for '{name}'.")
                time.sleep(sleep)
                pbar.update(1)

        return entrez_mapped

    @staticmethod
    def text2term_taxon_name2taxid(taxon_names: list[str], min_score=0.8) -> dict:
        """Map taxonomy names using text2term."""
        if not text2term.cache_exists("NCBITaxon"):
            print("Caching NCBITaxon ontology...")
            text2term.cache_ontology(
                ontology_url="http://purl.obolibrary.org/obo/ncbitaxon.owl",
                ontology_acronym="NCBITaxon",
            )

        taxon_names = list(set(taxon_names))

        with tqdm(desc="Text2term taxon mapping", unit="batch") as pbar:
            df_cached = text2term.map_terms(
                source_terms=taxon_names,
                target_ontology="NCBITaxon",
                use_cache=True,
                min_score=min_score,
                max_mappings=1,
            )
            pbar.update(1)

        text2term_mapped = dict(zip(df_cached["Source Term"], df_cached["Mapped Term CURIE"]))
        for name, taxid in text2term_mapped.items():
            text2term_mapped[name] = {
                "taxid": int(taxid.split(":")[1]),
                "mapping_tool": "text2term",
            }

        return text2term_mapped


class DiseaseServices:
    """Services for disease name to ID mapping."""

    @staticmethod
    def text2term_disease_name2id(
        disease_names: list[str],
        ontology: str = "MONDO",
        ontology_url: str = "http://purl.obolibrary.org/obo/mondo.owl",
    ) -> dict:
        """Map disease names using text2term."""
        if not text2term.cache_exists(ontology):
            print(f"Caching {ontology} ontology...")
            text2term.cache_ontology(ontology_url=ontology_url, ontology_acronym=ontology)

        with tqdm(desc=f"Text2term {ontology} mapping", unit="batch") as pbar:
            core_disease_map_df = text2term.map_terms(
                source_terms=list(set(disease_names)),
                target_ontology=ontology,
                use_cache=True,
                min_score=0.8,
                max_mappings=1,
            )
            pbar.update(1)

        if core_disease_map_df is None or core_disease_map_df.empty:
            return {}

        filtered_map_df = core_disease_map_df[
            ~core_disease_map_df["Mapped Term CURIE"]
            .astype(str)
            .str.contains("NCBITAXON", na=False)
        ]
        filtered_map_dict = dict(
            zip(filtered_map_df["Source Term"], filtered_map_df["Mapped Term CURIE"])
        )
        mapped = {
            name: {"id": _id, "mapping_tool": "text2term", "min_score": 0.8}
            for name, _id in filtered_map_dict.items()
        }
        return mapped


class ProteinServices:
    """Services for protein-related data retrieval."""

    @staticmethod
    async def get_protein_function(
        session: aiohttp.ClientSession, uniprot_id: str
    ) -> Optional[dict]:
        """Get protein function description from UniProt."""
        url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
        try:
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    for comment in data.get("comments", []):
                        if comment.get("commentType") == "FUNCTION":
                            texts = comment.get("texts", [])
                            if texts and texts[0].get("value"):
                                return {uniprot_id: {"description": texts[0]["value"]}}
                return None
        except Exception as e:
            tqdm.write(f"Error fetching {uniprot_id}: {str(e)}")
            return None

    @staticmethod
    async def get_batch_protein_functions(
        uniprot_ids: List[str], batch_size: int = 5, delay: float = 1.0
    ) -> List[dict]:
        """Get protein functions in batches with progress tracking."""
        results = []
        uniprot_ids = list(set(uniprot_ids))
        connector = aiohttp.TCPConnector(limit=batch_size)

        async with aiohttp.ClientSession(connector=connector) as session:
            with tqdm(
                desc="Fetching protein functions", total=len(uniprot_ids), unit="proteins"
            ) as pbar:
                for i in range(0, len(uniprot_ids), batch_size):
                    batch = uniprot_ids[i : i + batch_size]
                    tasks = [ProteinServices.get_protein_function(session, uid) for uid in batch]
                    batch_results = await asyncio.gather(*tasks)
                    results.extend([r for r in batch_results if r is not None])
                    pbar.update(len(batch))
                    await asyncio.sleep(delay)

        return results


class GeneServices:
    """Services for gene-related data retrieval."""

    @staticmethod
    async def get_gene_summary(session: aiohttp.ClientSession, gene_id: str) -> Optional[dict]:
        """Get gene summary from NCBI."""
        url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        params = {
            "db": "gene",
            "id": gene_id,
            "retmode": "json",
            "tool": "Microbiome Knowledge Graph",
            "email": "bazhang@scripps.edu",
            "api_key": os.getenv("NCBI_API_KEY"),
        }
        try:
            async with session.get(url, params=params) as resp:
                if resp.status == 200:
                    data = await resp.json()
                    summary = (
                        data.get("result", {})
                        .get(str(gene_id), {})
                        .get("summary", "No summary found.")
                    )
                    return {gene_id: {"description": summary}}
                return None
        except Exception as e:
            tqdm.write(f"Error fetching gene {gene_id}: {str(e)}")
            return None

    @staticmethod
    async def get_batch_gene_summaries(
        gene_ids: List[str], batch_size: int = 10, delay: float = 1.0
    ) -> List[dict]:
        """Get gene summaries in batches with progress tracking."""
        results = []
        connector = aiohttp.TCPConnector(limit=batch_size)

        async with aiohttp.ClientSession(connector=connector) as session:
            with tqdm(desc="Fetching gene summaries", total=len(gene_ids), unit="genes") as pbar:
                for i in range(0, len(gene_ids), batch_size):
                    batch = gene_ids[i : i + batch_size]
                    tasks = [GeneServices.get_gene_summary(session, gid) for gid in batch]
                    batch_results = await asyncio.gather(*tasks)
                    results.extend([r for r in batch_results if r is not None])
                    pbar.update(len(batch))
                    await asyncio.sleep(delay)

        return results


class GOServices:
    """Services for Gene Ontology data retrieval."""

    @staticmethod
    async def get_go_definitions(
        go_ids: List[str], batch_size: int = 200, delay: float = 0.25
    ) -> Dict[str, dict]:
        """Get GO term definitions with progress tracking."""
        BASE = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"

        def chunks(lst, n):
            for i in range(0, len(lst), n):
                yield lst[i : i + n]

        async def query(session: aiohttp.ClientSession, ids: List[str]) -> Dict[str, dict]:
            url = BASE + ",".join(ids)
            params = {"fields": "id,definition"}
            async with session.get(url, params=params, headers={"Accept": "application/json"}) as r:
                r.raise_for_status()
                js = await r.json()

            output_d = {}
            for term in js.get("results", []):
                go_id = term.get("id")
                descr = term.get("definition", {}).get("text")
                if go_id and descr:
                    output_d[go_id] = {"description": descr, "annotation_tool": "EBI_QuickGO"}
            return output_d

        results = {}
        timeout = aiohttp.ClientTimeout(total=None, connect=10, sock_read=30)

        go_chunks = list(chunks(go_ids, batch_size))
        async with aiohttp.ClientSession(timeout=timeout) as session:
            with tqdm(desc="Fetching GO definitions", total=len(go_chunks), unit="batches") as pbar:
                for group in go_chunks:
                    batch_results = await query(session, group)
                    results.update(batch_results)
                    pbar.update(1)
                    await asyncio.sleep(delay)

        return results


class UMLSClient:
    """UMLS API client for disease name to CUI mapping."""

    def __init__(self, api_key: str, max_concurrent: int = 10):
        self.api_key = api_key
        self.tgt_url: Optional[str] = None
        self.semaphore = asyncio.Semaphore(max_concurrent)

    async def get_tgt(self, session: aiohttp.ClientSession):
        """Get Ticket Granting Ticket."""
        data = {"apikey": self.api_key}
        async with session.post("https://utslogin.nlm.nih.gov/cas/v1/api-key", data=data) as resp:
            if resp.status != 201:
                raise RuntimeError(f"Failed to get TGT: {resp.status}")
            self.tgt_url = resp.headers["location"]

    async def get_st(self, session: aiohttp.ClientSession):
        """Get Service Ticket."""
        if not self.tgt_url:
            await self.get_tgt(session)
        data = {"service": "http://umlsks.nlm.nih.gov"}
        async with session.post(self.tgt_url, data=data) as resp:
            return await resp.text()

    async def get_cui(self, session: aiohttp.ClientSession, term: str) -> Optional[str]:
        """Get CUI for a term."""
        async with self.semaphore:
            try:
                st = await self.get_st(session)
                params = {"string": term, "ticket": st, "pageSize": 1, "searchType": "exact"}
                url = "https://uts-ws.nlm.nih.gov/rest/search/current"
                async with session.get(url, params=params) as resp:
                    if resp.status != 200:
                        return None
                    data = await resp.json()
                    results = data["result"]["results"]
                    return f"UMLS:{results[0]['ui']}" if results else ""
            except Exception as e:
                tqdm.write(f"Failed for {term}: {e}")
                return None

    async def query_cuis(self, terms: List[str]) -> Dict[str, dict]:
        """Query CUIs for multiple terms with progress tracking."""
        results = {}
        async with aiohttp.ClientSession() as session:
            if not self.tgt_url:
                await self.get_tgt(session)

            with tqdm(desc="UMLS CUI mapping", total=len(terms), unit="terms") as pbar:
                tasks = [self.get_cui(session, term) for term in terms]
                cuis = await asyncio.gather(*tasks)

                for term, cui in zip(terms, cuis):
                    results[term] = {"id": cui, "mapping_tool": "UMLS", "search_type": "exact"}
                    pbar.update(1)

        return results


class NCITServices:
    """Services for NCIT ontology data retrieval."""

    @staticmethod
    def get_ncit_taxon_description(taxon_names: list[str]) -> dict:
        """Get taxon descriptions from NCIT."""
        NCIT_API_KEY = os.getenv("NCIT_API_KEY")
        search_url = "https://data.bioontology.org/search"
        taxon_names = set(taxon_names)
        mapping_result = {}

        with tqdm(desc="NCIT taxon descriptions", total=len(taxon_names), unit="names") as pbar:
            for name in taxon_names:
                params = {
                    "q": name,
                    "ontologies": "NCIT",
                    "apikey": NCIT_API_KEY,
                }
                response = requests.get(search_url, params=params)
                data = response.json()

                for result in data.get("collection", []):
                    if result:
                        ncit_output = {
                            "name": result.get("prefLabel").lower(),
                            "description": f"{result.get('definition')[0]} [NCIT]"
                            if "definition" in result
                            else "",
                            "xrefs": {"ncit": result.get("@id").split("#")[1]},
                        }
                        if ncit_output["name"] == name:
                            mapping_result[name] = ncit_output
                            del mapping_result[name]["name"]
                pbar.update(1)

        return mapping_result


class BiothingsServices:
    """Services for biothings API data retrieval."""

    @staticmethod
    def get_taxon_info_from_bt(taxids: list) -> dict:
        """Get taxon information from biothings taxonomy API."""
        taxids = set(taxids)
        get_taxon = bt.get_client("taxon")

        with tqdm(desc="Fetching taxon info", unit="batch") as pbar:
            taxon_info = get_taxon.gettaxa(
                taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"]
            )
            pbar.update(1)

        taxon = {}
        for info in taxon_info:
            if "notfound" not in info.keys():
                taxon[info["query"]] = {
                    "id": f"NCBITaxon:{int(info['_id'])}",
                    "taxid": int(info["_id"]),
                    "name": info["scientific_name"],
                    "parent_taxid": int(info["parent_taxid"]),
                    "lineage": info["lineage"],
                    "rank": info["rank"],
                }
            else:
                taxon[info["query"]] = {
                    "id": f"NCBITaxon:{int(info['query'])}",
                    "taxid": int(info["query"]),
                }
        return taxon

    @staticmethod
    def bt_get_disease_info(ids: list) -> dict:
        """Get disease information from biothings disease API."""
        ids = set(ids)
        get_disease = bt.get_client("disease")

        with tqdm(desc="Fetching disease info", unit="batch") as pbar:
            d_queried = get_disease.querymany(
                ids,
                scopes=[
                    "mondo.mondo",
                    "mondo.xrefs.hp",
                    "mondo.xrefs.omim",
                    "mondo.xrefs.umls",
                    "mondo.xrefs.umls_cui",
                    "disease_ontology.xrefs.omim",
                    "disease_ontology.xrefs.umls_cui",
                ],
                fields=["mondo", "mondo.definition"],
            )
            pbar.update(1)

        d_info_all = {}
        for info in d_queried:
            query = info["query"]
            if "notfound" in info:
                continue

            current_score = info.get("_score", 0)
            existing = d_info_all.get(query)
            existing_score = existing.get("_score", -1) if existing else -1

            # determine prefix and id
            if current_score > existing_score:
                mondo_data = info.get("mondo", {})

                if query.isdigit():
                    prefix = "omim"
                    _id = f"OMIM:{query}"
                elif query.startswith("C") and query[1:].isdigit():
                    prefix = "umls"
                    _id = f"UMLS:{query}"
                else:
                    prefix = query.split(":")[0].lower() if ":" in query else "unknown"
                    _id = query

                d_info_all[_id] = {
                    "id": info["_id"],
                    "name": mondo_data.get("label", "").lower(),
                    "description": mondo_data.get("definition"),
                    "type": "biolink:Disease",
                    "xrefs": {prefix: _id if info["_id"] != _id else info["_id"]},
                    "_score": current_score,
                }

        # clean up
        for v in d_info_all.values():
            v.pop("_score", None)
            if v.get("description") is None:
                v.pop("description", None)

        return d_info_all


class PathwayServices:
    """Services for pathway data retrieval."""

    @staticmethod
    def get_smpdb_pathway_description() -> dict:
        """Get SMPDB pathway descriptions from CSV file."""

        from .reader import extract_file_from_zip

        zip_path = os.path.join("downloads", "smpdb_pathways.csv.zip")
        smpdb_csv = extract_file_from_zip(zip_path, "smpdb_pathways.csv")

        with tqdm(desc="Loading SMPDB pathways", unit="file") as pbar:
            smpdb_df = pd.read_csv(smpdb_csv, usecols=["SMPDB ID", "Name", "Description"])
            pbar.update(1)

        return {
            name.lower(): {
                "smpdb": f"SMPDB:{sid}",
                "description": descr,
            }
            for sid, name, descr in zip(
                smpdb_df["SMPDB ID"], smpdb_df["Name"], smpdb_df["Description"]
            )
        }


class UberonService:
    """Handles anatomy EBI OLS service."""

    async def async_query_anatomical_entity_to_uberon_id(
        self, session, term, url, match_type="exact", ontology="uberon", rows=1
    ):
        """Async helper function that now accepts a URL."""
        params = {"q": term, "ontology": ontology, "rows": rows, "start": 0}

        if match_type == "exact":
            params["exact"] = "true"
        elif match_type == "partial":
            params["exact"] = "false"
            params["queryFields"] = "label"
        elif match_type == "fuzzy":
            params["q"] = term + "~"
            params["exact"] = "false"
        else:
            raise ValueError("match_type must be 'exact', 'partial', or 'fuzzy'.")

        try:
            async with session.get(url, params=params) as resp:
                resp.raise_for_status()
                data = await resp.json()

                for doc in data.get("response", {}).get("docs", []):
                    iri = doc.get("iri")
                    if iri and "UBERON_" in iri:
                        uberon_id = iri.split("/")[-1].replace("_", ":")

                        return term, {
                            "id": uberon_id,
                            "name": doc.get("label"),
                            "original_name": term,
                            "type": "biolink:AnatomicalEntity",
                        }
        except aiohttp.ClientError as e:
            print(f"❗️Error fetching term '{term}': {e}")
        return term, None

    async def async_anatomical_entities_to_uberon_ids(
        self, terms, match_type="exact", base_url="https://www.ebi.ac.uk/ols4/api/search"
    ):
        """
        Query OLS for a list of terms, using the provided base_url.
        """
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.async_query_anatomical_entity_to_uberon_id(
                    session, term, url=base_url, match_type=match_type
                )
                for term in terms
            ]

            all_results = await tqdm.gather(*tasks, desc="Querying UBERON IDs...")

            return {term: info for term, info in all_results if info is not None}

    def async_run_anatomical_entities_to_uberon_ids(self, terms, match_type="exact"):
        base_url = "https://www.ebi.ac.uk/ols4/api/search"
        results = asyncio.run(
            self.async_anatomical_entities_to_uberon_ids(
                terms, match_type=match_type, base_url=base_url
            )
        )
        print("🎉 UBERON ID mapping completed!")
        return results
