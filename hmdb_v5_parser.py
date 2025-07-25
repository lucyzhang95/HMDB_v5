import asyncio
import json
import os
import pathlib
import pickle
import re
import ssl
import time
import uuid
import zipfile
from itertools import chain
from typing import Dict, Iterator, List, Optional, Union

import aiohttp
import biothings_client as bt
import pandas as pd
import requests
import text2term
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from lxml import etree as ET
from tqdm.auto import tqdm

load_dotenv()
CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    """
    :param obj:
    :param f_name: files should only be existing in the cache directory
    :return:
    """
    with open(os.path.join(CACHE_DIR, f_name), "wb") as out_f:
        pickle.dump(obj, out_f)


def load_pickle(f_name):
    """

    :param f_name:
    :return:
    """
    path = os.path.join(CACHE_DIR, f_name)
    return pickle.load(open(path, "rb")) if os.path.exists(path) else None


def save_json(obj, f_name):
    with open(os.path.join(CACHE_DIR, f_name), "w") as out_f:
        json.dump(obj, out_f, indent=4)


def extract_file_from_zip(zip_path, expected_filename):
    extract_dir = os.path.dirname(zip_path)
    extracted_path = os.path.join(extract_dir, expected_filename)

    if not os.path.isfile(extracted_path):
        with zipfile.ZipFile(zip_path, "r") as zip_f:
            if expected_filename in zip_f.namelist():
                zip_f.extract(expected_filename, path=extract_dir)
            else:
                raise FileNotFoundError(f"{expected_filename} not found in {zip_path}")

    return extracted_path


def strip_tag_namespace(tag: str) -> str:
    """Strip namespace from main element tag
    Each element tag from HMDB contains "http://www.hmdb.ca" as namespace
    e.g. {http://www.hmdb.ca}accession
    Remove the namespace for easier data processing

    :param tag: element tags (e.g. <accession>HMDB0000001</accession>)
    :return: original tag without namespace
    """
    return tag.split("}", 1)[-1] if "}" in tag else tag


def get_all_microbe_names(input_xml: str | pathlib.Path) -> Iterator[str]:
    """Extracts microbe taxon names
    metabolite > ontology > root (2nd) > Disposition > descendants > descendant > Microbe
    <term>Microbe</term> is in the 2nd root under ontology

    :param input_xml: "downloads/hmdb_metabolites.xml"
    :return: yields all microbe taxon names
    """
    ns_uri = "http://www.hmdb.ca"
    term_tag = f"{{{ns_uri}}}term"
    descendant_tag = f"{{{ns_uri}}}descendant"
    ontology_tag = f"{{{ns_uri}}}ontology"

    for _, elem in ET.iterparse(input_xml, events=("end",), tag=f"{{{ns_uri}}}metabolite"):
        ontology = elem.find(ontology_tag)
        if ontology is not None:
            for root in ontology.findall(f"{f'{{{ns_uri}}}root'}"):
                term_elem = root.find(term_tag)
                if term_elem is not None and term_elem.text == "Disposition":
                    for descendant in root.findall(f".//{descendant_tag}"):
                        term = descendant.find(term_tag)
                        if term is not None and term.text == "Microbe":
                            microbe_terms = descendant.findall(f".//{term_tag}")
                            for mt in microbe_terms[1:]:  # skip first term == Microbe
                                if mt.text:
                                    yield mt.text.strip().lower()
        elem.clear()


def get_all_diseases(input_xml: Union[str, pathlib.Path]) -> dict[str, str]:
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    disease2omim = {}
    for metabolite in root.findall("hmdb:metabolite", namespace):
        diseases_elem = metabolite.find("hmdb:diseases", namespace)
        if diseases_elem is None:
            continue

        for disease_elem in diseases_elem.findall("hmdb:disease", namespace):
            name_elem = disease_elem.find("hmdb:name", namespace)
            omim_elem = disease_elem.find("hmdb:omim_id", namespace)
            if name_elem is not None and name_elem.text:
                disease_name = name_elem.text.strip().lower()
                omim_id = (
                    omim_elem.text.strip() if omim_elem is not None and omim_elem.text else None
                )
                disease2omim[disease_name] = int(omim_id) if omim_id else None

    return disease2omim


def ete3_taxon_name2taxid(taxon_names: list) -> dict:
    """Use ete3 to map taxonomy names to NCBI taxonomy ids
    ete3 is good at mapping exact taxonomy names and fast without accessing the API

    :param taxon_names: a list of taxon names
    :return: a dictionary mapping taxon names to taxids
    e.g., {'veillonella alcalescens': {'taxid': 29466, 'mapping_tool': 'ete3'}, ...}
    """
    ete3_mapped = {}
    taxon_names = set(taxon_names)

    if not NCBITaxa:
        ssl._create_default_https_context = ssl._create_unverified_context
        ncbi = NCBITaxa()
        ncbi.update_taxonomy_database()

    ncbi = NCBITaxa()
    ete3_name2taxid = ncbi.get_name_translator(taxon_names)
    for name, taxid in ete3_name2taxid.items():
        if taxid:
            ete3_mapped[name] = {"taxid": int(taxid[0]), "mapping_tool": "ete3"}
    return ete3_mapped


def entrez_taxon_name2taxid(
    taxon_names: list[str],
    email: str,
    sleep: float = 0.34,
    retries: int = 3,
    backoff_factor: int = 2,
) -> dict:
    """Map taxonomy names to NCBI taxonomy ids using the Entrez API with retry logic.

    :param taxon_names: A list of taxonomy names to query.
    :param email: Your email address for NCBI Entrez.
    :param sleep: The base time to sleep between different queries.
    :param retries: The maximum number of retry attempts for a failed query.
    :param backoff_factor: The factor by which to increase the delay between retries.
    :return: A dictionary mapping taxon names to taxids.
    """
    Entrez.email = email
    entrez_mapped = {}

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
                print(f"Entrez query failed for '{name}' on attempt {attempt + 1}/{retries}: {e}")
                if attempt < retries - 1:
                    print(f"Retrying in {delay} seconds...")
                    time.sleep(delay)
                    delay *= backoff_factor
                else:
                    print(f"All retry attempts failed for '{name}'.")
        time.sleep(sleep)

    return entrez_mapped


def text2term_taxon_name2taxid(taxon_names: list[str], min_score=0.8) -> dict:
    """

    :param taxon_names:
    :param min_score:
    :return:
    e.g., {'pseudomonas uorescens': {'taxid': 294, 'mapping_tool': 'text2term'}, ...}
    """
    if not text2term.cache_exists("NCBITaxon"):
        text2term.cache_ontology(
            ontology_url="http://purl.obolibrary.org/obo/ncbitaxon.owl",
            ontology_acronym="NCBITaxon",
        )

    taxon_names = list(set(taxon_names))
    df_cached = text2term.map_terms(
        source_terms=taxon_names,
        target_ontology="NCBITaxon",
        use_cache=True,
        min_score=min_score,
        max_mappings=1,
    )

    text2term_mapped = dict(zip(df_cached["Source Term"], df_cached["Mapped Term CURIE"]))
    for name, taxid in text2term_mapped.items():
        text2term_mapped[name] = {"taxid": int(taxid.split(":")[1]), "mapping_tool": "text2term"}
    return text2term_mapped


def manual_correct_text2term_map(text2term_mapped: dict) -> dict:
    manual_corrections = {
        "geobacillus thermoglucosidasius": 1426,
        "bacteroides spp.": 816,
        "clostridium difficile": 1496,
        "pseudomonas pseudomaleii": 28450,
        "lactobacillus plantarum": 1590,
        "clostridium butylicum": 1492,
        "corynebacterium jekeium": 38289,
    }

    for name, corrected_taxid in manual_corrections.items():
        if name in text2term_mapped:
            text2term_mapped[name]["taxid"] = corrected_taxid
            text2term_mapped[name]["mapping_tool"] = "manual"

    del text2term_mapped["gram-negative bacteria"]
    return text2term_mapped


def manual_taxon_name2taxid(taxon_names: list[str]) -> dict:
    """

    :param taxon_names:
    :return:
    e.g., {'clostridia propionicum': {'taxid': 28446, 'mapping_tool': 'manual'}, ...}
    """
    raw_mapping = {
        "rlzodopseudomonas spheroides": 1063,
        "clostridium calortolerans": 36832,
        "clostridium felsenium": 36839,
        "rhodobacter spaeroides": 1063,
        "clostridium propylbutyricum": 1485,
        "algibacter sp. aqp096": 1872428,
        "pseudomonas ligustri": 317,
        "mycobacterium smegmatis": 1772,
        "clostridium stricklandii": 1511,
        "clostridium species": 1485,
        "biﬁdobacterium": 1678,
        "pseudomonas sp. dsm 2874": 306,
        "methanothrix sochngenii": 2223,
        "pseudomonas sp. st-200": 306,
        "clostridium lituseburense": 1537,
        "muricauda lutaonensis": 516051,
        "clostridia propionicum": 28446,
        "citrobacter frundii": 546,
        "meningococcus": 487,
        "pseudomonas orvilla": 303,
        "clostridium xiva": 543317,
        "chromobacterium prodigiosum": 615,
        "clostridium mangenoyi": 1540,
        "clostridium glycolycum": 36841,
        "streptococcus group b": 1319,
    }

    manual_mapped = {
        name: {"taxid": raw_mapping[name], "mapping_tool": "manual"}
        for name in list(set(taxon_names))
        if name in raw_mapping
    }
    return manual_mapped


def get_taxon_info_from_bt(taxids) -> dict:
    """

    :param taxids:
    :return:
        {
       "28450":{
          "id":"NCBITaxon:28450",
          "taxid":28450,
          "name":"burkholderia pseudomallei",
          "parent_taxid":111527,
          "lineage":[
             28450,
             111527,
             32008,
             119060,
             80840,
             28216,
             1224,
             3379134,
             2,
             131567,
             1
          ],
          "rank":"species"
       }"..."
    }
    """
    taxids = set(taxids)
    get_taxon = bt.get_client("taxon")
    taxon_info = get_taxon.gettaxa(
        taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"]
    )

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


def get_ncit_taxon_description(taxon_names):
    """

    :param taxon_names:
    :return:
    {'serratia': {'description':
    'A genus of small motile peritrichous bacteria in the Enterobacteriacaea family
    consisting of Gram-negative rods.
    [NCIT]',
    'xrefs': {'ncit': 'C86010', }} ...}
    """
    NCIT_API_KEY = os.getenv("NCIT_API_KEY")
    search_url = "https://data.bioontology.org/search"
    taxon_names = set(taxon_names)
    mapping_result = {}
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
    return mapping_result


def add_description2taxon_info(taxon_info: dict, descriptions: dict) -> dict:
    for _, info in taxon_info.items():
        name = info.get("name")
        if name in descriptions:
            info.update(descriptions[name])
        else:
            pass

    return taxon_info


def get_full_taxon_info(mapped_taxon_names: dict, taxon_info: dict) -> dict:
    """A dictionary of taxon metadata keyed by taxon name,
    using taxid lookup from `taxon_info`.

    :param mapped_taxon_names: Mapping of taxon name to a dict with taxon taxid, lineage, rank...
    {'clostridia propionicum': {'taxid': 28446, 'mapping_tool': 'manual'} ...}
    :param taxon_info: Mapping of taxid (as str) to a dictionary with mapped taxon names and taxid.
    '28450': {'id': 'NCBITaxon:28450', 'taxid': 28450, 'name': 'burkholderia pseudomallei', 'parent_taxid': 111527,
    'lineage': [28450, 111527, 32008, 119060, 80840, 28216, 1224, 3379134, 2, 131567, 1], 'rank': 'species',
    'description': 'A species of aerobic, Gram-negative, rod shaped bacteria assigned to the phylum Proteobacteria.
    This species is motile, non-spore forming, oxidase and catalase positive and indole negative.
    B. pseudomallei is found in contaminated soil, water,
    and produce and causes melioidosis in humans with the highest rate of disease occurring in southeast Asia.
    [NCIT]', 'xrefs': {'ncit': 'C86010' }} ...}

    :return:

    """
    full_taxon = {
        name: taxon_info[str(taxon_entry["taxid"])]
        for name, taxon_entry in mapped_taxon_names.items()
        if str(taxon_entry["taxid"]) in taxon_info
    }
    return full_taxon


def get_cuis_sync(api_key: str, disease_names: List[str]):
    client = UMLSClient(api_key=api_key)

    async def inner():
        return await client.query_cuis(disease_names)

    return asyncio.run(inner())


def text2term_disease_name2id(
    disease_names, ontology="MONDO", ontology_url="http://purl.obolibrary.org/obo/mondo.owl"
):
    """

    :param disease_names:
    :param ontology:
    :param ontology_url:
    {"EFO": "http://www.ebi.ac.uk/efo/efo.owl", "NCIT": "http://purl.obolibrary.org/obo/ncit.owl"}
    :return:
    """
    if not text2term.cache_exists(ontology):
        text2term.cache_ontology(ontology_url=ontology_url, ontology_acronym=ontology)

    core_disease_map_df = text2term.map_terms(
        source_terms=list(set(disease_names)),
        target_ontology=ontology,
        use_cache=True,
        min_score=0.8,
        max_mappings=1,
    )

    if core_disease_map_df is None or core_disease_map_df.empty:
        return {}

    filtered_map_df = core_disease_map_df[
        ~core_disease_map_df["Mapped Term CURIE"].astype(str).str.contains("NCBITAXON", na=False)
    ]
    filtered_map_dict = dict(
        zip(filtered_map_df["Source Term"], filtered_map_df["Mapped Term CURIE"])
    )
    mapped = {
        name: {"id": _id, "mapping_tool": "text2term", "min_score": 0.8}
        for name, _id in filtered_map_dict.items()
    }
    return mapped


def manual_disease_name2id(disease_names: list[str]) -> dict:
    raw_mapping = {
        "early preeclampsia": "MONDO:0005081",  # preeclampsia
        "late-onset preeclampsia": "MONDO:0005081",  # preeclampsia
        "perillyl alcohol administration for cancer treatment": "GO:0018457",  # perillyl-alcohol dehydrogenase (NAD+) activity
        "3-hydroxyisobutyric acid dehydrogenase deficiency": "MONDO:0009371",  # 3-hydroxyisobutyric aciduria
        "attachment loss": "UMLS:C0206114",  # periodontal attachment loss
        "periodontal probing depth": "UMLS:C1882338",  # periodontal probing
        "methylmalonic aciduria mitochondrial encephelopathy leigh-like": "UMLS:C1855119",  # methylmalonic aciduria
        "functional hypothalamic amenorrhea": "UMLS:C0341862",  # hypothalamic amenorrhea
        "peroxisomal disorders, new type, liver": "UMLS:C5568675",  # liver disease due to peroxisomal disease
        "dopamine-serotonin vesicular transport defect": "MONDO:0018130",  # brain dopamine-serotonin vesicular transport disease
        "prosthesis/missing teeth": "UMLS:C0080233",  # tooth Loss
        "small intestinal malabsorption": "UMLS:C1833057",  # malabsorption (small intestine)
        "refractory localization-related epilepsy": "UMLS:C0472349",  # localization-related symptomatic epilepsy
        "d-glyceric acidura": "MONDO:0009070",  # d-glyceric aciduria
        "methyl formate exposure": "ECTO:9001470",  # exposure to methyl formate
        "formic acid intoxication": "ECTO:9000376",  # exposure to formic acid
        "idiopathic oro-facial pain": "MONDO:0018362",  # persistent idiopathic facial pain
        "serine deficiency syndrome, infantile": "MONDO:0035004",  # serine biosynthesis pathway deficiency, infantile/juvenile form
        "hepatic and biliary malignancies": "MONDO:0002514",  # hepatobiliary neoplasm
        "acute seizures": "UMLS:C0036572",  # seizures
        "methamphetamine (map) psychosis": "MONDO:0005465",  # methamphetamine-induced psychosis
        "nicotinamide adenine dinucleotide deficiency": "UMLS:C1283629",  # deficiency of NAD+ nucleosidase
        "homozygous sickle cell disease": "MONDO:0011382",  # sickle cell disease
        "prepartum depression": "UMLS:C0011570",  # mental depressionm -> only postpartum depression exists
        "nucleotide depletion syndrome": "MONDO:0018158",  # mitochondrial DNA depletion syndrome
        "terminal aldosterone biosynthesis defects": "MONDO:0018541",  # familial hypoaldosteronism (Aldosterone synthase deficiency is a rare inherited defect of the final step of aldosterone biosynthesis)
        "glutaryl-coa dehydrogenase deficiency (gdhd)": "MONDO:0009281",  # glutaryl-CoA dehydrogenase deficiency
        "cancer with metastatic bone disease": "UMLS:C5444038",  # metastatic bone disease
        "neuroinfection": "UMLS:C0870953",  # neuroinfections
        "dermal fibroproliferative disorder": "UMLS:C1304434",  # dermal elastolysis
        "d-lactic acidosis and short bowel syndrome": "MONDO:0015183",  # short bowel syndrome
        "tert-amyl-methyl ether exposed": "ECTO:9000278",  # exposure to ether
        "vessel occlusion": "MONDO:0020673",  # arterial occlusion
        "cresol poisoning ibs": "ECTO:9001047",  # exposure to cresol
        "dimethyl sulfide poisoning": "UMLS:C2062726",  # poisoning by sulfides
        "quetiapine poisoning": "ECTO:9000356",  # exposure to quetiapine
    }

    manual_mapped = {
        name: {"id": raw_mapping[name], "mapping_tool": "manual"}
        for name in list(set(disease_names))
        if name in raw_mapping
    }
    return manual_mapped


def bt_get_disease_info(ids):
    ids = set(ids)
    get_disease = bt.get_client("disease")
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

            if re.fullmatch(r"\d+", query):
                prefix = "omim"
                _id = f"OMIM:{query}"
            elif re.fullmatch(r"C\d+", query):
                prefix = "umls"
                _id = f"UMLS:{query}"
            else:
                prefix = query.split(":")[0].lower()
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


async def get_protein_function(session, uniprot_id):
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
                print({uniprot_id: "Function not found."})
            else:
                print(f"Failed: HTTP {response.status}")
    except Exception as e:
        print(uniprot_id, f"Error: {str(e)}")


async def get_batch_protein_functions(uniprot_ids: List[str], batch_size=5, delay=1.0):
    results = []
    uniprot_ids = list(set(uniprot_ids))
    connector = aiohttp.TCPConnector(limit=batch_size)
    async with aiohttp.ClientSession(connector=connector) as session:
        for i in range(0, len(uniprot_ids), batch_size):
            batch = uniprot_ids[i : i + batch_size]
            tasks = [get_protein_function(session, uid) for uid in batch]
            batch_results = await asyncio.gather(*tasks)
            results.extend([r for r in batch_results if r is not None])  # filter out None results
            await asyncio.sleep(delay)
    return results


def get_all_uniprot_ids_from_hmdb(input_xml) -> dict[str, str]:
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    protein2uniport = {}
    for metabolite in root.findall("hmdb:metabolite", namespace):
        protein_assoc_elem = metabolite.find("hmdb:protein_associations", namespace)
        if protein_assoc_elem is None:
            continue

        for protein_elem in protein_assoc_elem.findall("hmdb:protein", namespace):
            name_elem = protein_elem.find("hmdb:name", namespace)
            uniprot_elem = protein_elem.find("hmdb:uniprot_id", namespace)
            if name_elem is not None and name_elem.text:
                protein_name = name_elem.text.strip().lower()
                uniprot_id = (
                    uniprot_elem.text.strip()
                    if uniprot_elem is not None and uniprot_elem.text
                    else None
                )
                protein2uniport[protein_name] = uniprot_id if uniprot_id else None
    return protein2uniport


def get_all_uniprot_ids_from_hmdbp(input_xml) -> list:
    """Extracts all UniProt IDs from HMDBP XML file."""
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    uniprot_ids = [
        elem.text.strip()
        for protein in root.findall("hmdb:protein", namespace)
        if (elem := protein.find("hmdb:uniprot_id", namespace)) is not None
        and elem.text
        and elem.text.strip()
    ]
    return uniprot_ids


def uniprot_id2entrezgene(uniprot_ids: list[str]) -> dict:
    """

    :param uniprot_ids:
    :return:
    e.g., {'Q8IYK8': {'gene_id': 'NCBIGene:161253', 'mapping_tool': 'bt'},...}
    """
    uniprot_ids = list(set(uniprot_ids))
    get_gene = bt.get_client("gene")
    gene_q = get_gene.querymany(uniprot_ids, scopes=["uniprot", "uniprot.Swiss-Prot"])

    entrezgene_mapped = {}
    for info in gene_q:
        if "notfound" in info:
            continue
        if "entrezgene" in info:
            entrezgene_mapped[info["query"]] = {
                "gene_id": f"NCBIGene:{info['entrezgene']}",
                "mapping_tool": "bt",
            }

    return entrezgene_mapped


async def get_gene_summary(session, gene_id):
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
                    data.get("result", {}).get(str(gene_id), {}).get("summary", "No summary found.")
                )
                return {gene_id: {"description": summary}}
            print({gene_id, f"Failed: HTTP {resp.status}"})
    except Exception as e:
        print(gene_id, f"Error: {str(e)}")


async def get_batch_gene_summaries(gene_ids: List[str], batch_size=10, delay=1.0):
    results = []
    connector = aiohttp.TCPConnector(limit=batch_size)
    async with aiohttp.ClientSession(connector=connector) as session:
        for i in range(0, len(gene_ids), batch_size):
            batch = gene_ids[i : i + batch_size]
            tasks = [get_gene_summary(session, gid) for gid in batch]
            batch_results = await asyncio.gather(*tasks)
            results.extend(batch_results)
            await asyncio.sleep(delay)
    return results


def get_smpdb_pathway_description():
    zip_path = os.path.join("downloads", "smpdb_pathways.csv.zip")
    smpdb_csv = extract_file_from_zip(zip_path, "smpdb_pathways.csv")
    smpdb_df = pd.read_csv(smpdb_csv, usecols=["SMPDB ID", "Name", "Description"])
    return {
        name.lower(): {
            "smpdb": f"SMPDB:{sid}",
            "description": descr,
        }
        for sid, name, descr in zip(smpdb_df["SMPDB ID"], smpdb_df["Name"], smpdb_df["Description"])
    }


async def get_go_definitions(
    go_ids: List[str], batch_size: int = 200, delay: float = 0.25
) -> Dict[str, str]:
    """

    :param go_ids:
    :param batch_size:
    :param delay:
    :return:
    """
    BASE = "https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/"

    def chunks(lst, n):
        for i in range(0, len(lst), n):
            yield lst[i : i + n]

    async def query(session: aiohttp.ClientSession, ids: List[str]) -> Dict[str, str]:
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

    async with aiohttp.ClientSession(timeout=timeout) as session:
        for group in chunks(go_ids, batch_size):
            results.update(await query(session, group))
            await asyncio.sleep(delay)

    return results


def get_all_go_terms_from_hmdbp(input_xml):
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    go_terms = set()
    for protein in root.findall("hmdb:protein", namespace):
        go_elem = protein.find("hmdb:go_classifications", namespace)
        if go_elem is None:
            continue
        for go_class in go_elem.findall("hmdb:go_class", namespace):
            go_id_elem = go_class.find("hmdb:go_id", namespace)
            if go_id_elem is not None and go_id_elem.text:
                go_terms.add(go_id_elem.text.strip())
    return list(go_terms)


def get_organism_type(node) -> str:
    """
    Inspect node['lineage'] for known taxids.
    Return the matching biolink CURIE, or Other if no match.
    Types include: 3 domains of life (Bacteria, Archaea, Eukaryota) and Virus.
    """
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


# TODO: Need to write a function check if cache files exist
def cache_metabolite_data(input_xml):
    # cache mapped taxon
    microbe_names = get_all_microbe_names(input_xml)
    microbes4query = list(set(microbe_names))
    save_pickle(microbes4query, "hmdb_v5_microbe_names.pkl")
    ete3_mapped = ete3_taxon_name2taxid(microbes4query)
    save_pickle(ete3_mapped, "ete3_name2taxid.pkl")

    no_hits = [name for name in microbes4query if name not in ete3_mapped]
    entrez_mapped = entrez_taxon_name2taxid(no_hits, email=os.getenv("EMAIL_ADDRESS"))
    save_pickle(entrez_mapped, "entrez_name2taxid.pkl")

    no_hits2 = [
        name for name in microbes4query if name not in ete3_mapped and name not in entrez_mapped
    ]
    text2term_mapped = text2term_taxon_name2taxid(no_hits2)
    text2term_mapped = manual_correct_text2term_map(text2term_mapped)
    save_pickle(text2term_mapped, "text2term_name2taxid.pkl")

    no_hits3 = [
        name
        for name in set(microbe_names)
        if name not in ete3_mapped and name not in entrez_mapped and name not in text2term_mapped
    ]
    manual_mapped = manual_taxon_name2taxid(no_hits3)
    save_pickle(manual_mapped, "manual_name2taxid.pkl")

    all_mapped_taxon_names = ete3_mapped | entrez_mapped | text2term_mapped | manual_mapped
    save_pickle(all_mapped_taxon_names, "all_taxon_name2taxid.pkl")

    all_mapped_taxon_cached = load_pickle("all_taxon_name2taxid.pkl")
    taxid2taxon = [int(taxid["taxid"]) for name, taxid in all_mapped_taxon_cached.items()]
    taxon_info = get_taxon_info_from_bt(taxid2taxon)
    taxon_sci_names = [info["name"] for info in taxon_info.values() if "name" in info]
    taxon_descr = get_ncit_taxon_description(taxon_sci_names)
    taxon_info_descr = add_description2taxon_info(taxon_info, taxon_descr)
    full_taxon_info = get_full_taxon_info(all_mapped_taxon_cached, taxon_info_descr)
    save_pickle(full_taxon_info, "original_taxon_name2taxid.pkl")

    # cache mapped diseases
    diseases = get_all_diseases(input_xml)
    save_pickle(diseases, "hmdb_v5_diseases.pkl")
    disease4text2term = [di_name for di_name, omim in diseases.items() if omim is None]
    text2term_mapped = text2term_disease_name2id(disease4text2term)
    save_pickle(text2term_mapped, "text2term_disease_name2id.pkl")
    disease2cui = [di_name for di_name in disease4text2term if di_name not in text2term_mapped]
    umls_mapped = get_cuis_sync(os.getenv("UMLS_API_KEY"), disease2cui)
    save_pickle(umls_mapped, "umls_disease_name2id.pkl")
    di_no_hit = [name for name, mapped in umls_mapped.items() if not mapped.get("id")]
    di_manual_mapped = manual_disease_name2id(di_no_hit)
    save_pickle(di_manual_mapped, "manual_disease_name2id.pkl")
    hmdb_di_mapped = {
        di_name: {"id": f"OMIM:{omim}", "mapping_tool": "hmdb_v5"}
        for di_name, omim in diseases.items()
        if omim
    }
    all_mapped_disease_names = hmdb_di_mapped | text2term_mapped | umls_mapped | di_manual_mapped
    filtered_all_mapped_disease_names = {
        name: mapped for name, mapped in all_mapped_disease_names.items() if not mapped["id"] == ""
    }
    save_pickle(filtered_all_mapped_disease_names, "all_disease_name2id.pkl")

    # TODO: need to check if the file exists before loading
    all_di_mapped = load_pickle("all_disease_name2id.pkl")

    di_ids = [
        mapped["id"]
        if "UMLS" not in mapped["id"] and "OMIM" not in mapped["id"]
        else mapped["id"].split(":")[1]
        for _, mapped in all_di_mapped.items()
    ]

    bt_di_info = bt_get_disease_info(di_ids)
    di_info = {}
    for name, mapped in all_di_mapped.items():
        _id = mapped["id"]
        if _id in bt_di_info:
            info = bt_di_info[_id].copy()
            info["original_name"] = name
            di_info[name] = info
    other_di = {
        name: {"id": mapped["id"], "name": name, "type": "biolink:Disease"}
        for name, mapped in all_di_mapped.items()
        if mapped["id"] not in bt_di_info
    }

    original_di_name_all = di_info | other_di
    save_pickle(original_di_name_all, "original_disease_name2id.pkl")

    # cache mapped proteins and protein functions
    mapped_proteins = get_all_uniprot_ids_from_hmdb(input_xml)
    save_pickle(mapped_proteins, "all_protein_name2uniprot.pkl")
    uniprot_ids = [uniprot for name, uniprot in mapped_proteins.items()]
    mapped_protein_descr = asyncio.run(get_batch_protein_functions(uniprot_ids))
    save_pickle(mapped_protein_descr, "uniprot_protein_functions.pkl")

    # cache entrezgene ids and gene summaries
    uniprot2entrez = uniprot_id2entrezgene(uniprot_ids)
    save_pickle(uniprot2entrez, "bt_uniprot2entrezgene.pkl")
    entrezgenes = [
        mapped["gene_id"].split(":")[1]
        for _, mapped in uniprot2entrez.items()
        if "gene_id" in mapped
    ]
    gene_descr = asyncio.run(get_batch_gene_summaries(entrezgenes))
    save_pickle(gene_descr, "entrezgene_summaries.pkl")

    # cache SMPDB pathway descriptions
    smpdb_pathway_descr = get_smpdb_pathway_description()
    save_pickle(smpdb_pathway_descr, "smpdb_pathway_descriptions.pkl")


def cache_protein_data(input_xml):
    # cache HMDBP protein functions
    hmdbp_uniprot_ids = get_all_uniprot_ids_from_hmdbp(input_xml)
    hmdbp_prot_func = asyncio.run(get_batch_protein_functions(hmdbp_uniprot_ids))
    save_pickle(hmdbp_prot_func, "hmdbp_uniprot_protein_functions.pkl")

    # cache gene summaries from HMDBP uniprot ids
    entrezgene2uniprot = uniprot_id2entrezgene(hmdbp_uniprot_ids)
    save_pickle(entrezgene2uniprot, "hmdbp_uniprot2entrezgene.pkl")
    entrezgenes = [
        gene_info["gene_id"].split(":")[1] for uniprot, gene_info in entrezgene2uniprot.items()
    ]
    hmdbp_gene_descr = asyncio.run(get_batch_gene_summaries(entrezgenes))
    save_pickle(hmdbp_gene_descr, "hmdbp_entrezgene_summaries.pkl")

    # cache GO terms from HMDBP
    go_terms = get_all_go_terms_from_hmdbp(input_xml)
    go_descr = asyncio.run(get_go_definitions(go_terms))
    save_pickle(go_descr, "hmdbp_go_definitions.pkl")


class UMLSClient:
    def __init__(self, api_key: str, max_concurrent: int = 10):
        self.api_key = api_key
        self.tgt_url: Optional[str] = None
        self.semaphore = asyncio.Semaphore(max_concurrent)

    async def get_tgt(self, session: aiohttp.ClientSession):
        data = {"apikey": self.api_key}
        async with session.post("https://utslogin.nlm.nih.gov/cas/v1/api-key", data=data) as resp:
            if resp.status != 201:
                raise RuntimeError(f"Failed to get TGT: {resp.status}")
            self.tgt_url = resp.headers["location"]

    async def get_st(self, session: aiohttp.ClientSession):
        if not self.tgt_url:
            await self.get_tgt(session)
        data = {"service": "http://umlsks.nlm.nih.gov"}
        async with session.post(self.tgt_url, data=data) as resp:
            return await resp.text()

    async def get_cui(self, session: aiohttp.ClientSession, term: str) -> Optional[str]:
        async with self.semaphore:
            try:
                st = await self.get_st(session)
                params = {"string": term, "ticket": st, "pageSize": 1, "searchType": "exact"}
                url = "https://uts-ws.nlm.nih.gov/rest/search/current"
                async with session.get(url, params=params) as resp:
                    if resp.status != 200:
                        print(f"Error querying {term}: status {resp.status}")
                        return None
                    data = await resp.json()
                    results = data["result"]["results"]
                    return f"UMLS:{results[0]['ui']}" if results else ""
            except Exception as e:
                print(f"Failed for {term}: {e}")
                return None

    async def query_cuis(self, terms: List[str]) -> Dict[str, Optional[str]]:
        results = {}
        async with aiohttp.ClientSession() as session:
            if not self.tgt_url:
                await self.get_tgt(session)
            tasks = [self.get_cui(session, term) for term in terms]
            cuis = await asyncio.gather(*tasks)
            results = {
                term: {"id": cui, "mapping_tool": "UMLS", "search_type": "exact"}
                for term, cui in zip(terms, cuis)
            }
        return results


class XMLParseHelper:
    namespace = {"hmdb": "http://www.hmdb.ca"}

    def __init__(self, input_xml):
        self.input_xml = input_xml

    def get_text(self, elem, tag):
        child = elem.find(f"hmdb:{tag}", self.namespace)
        return child.text.strip() if child is not None and child.text else None

    def get_list(self, elem, tag):
        return [e.text.lower() for e in elem.findall(f"hmdb:{tag}", self.namespace) if e.text]

    def remove_empty_none_values(self, obj):
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


class HMDB_Metabolite_Parse(XMLParseHelper):
    def __init__(self, input_xml):
        super().__init__(input_xml)
        self.parenthetical_pattern = re.compile(r"([^.?!]*?)\s*\(([^)]*?)\)")
        self.cached_taxon_info = load_pickle("original_taxon_name2taxid.pkl")
        self.cached_disease_info = load_pickle("original_disease_name2id.pkl")
        self.cached_protein_function = load_pickle("uniprot_protein_functions.pkl")
        self.cached_pathway_descr = load_pickle("smpdb_pathway_descriptions.pkl")

    def get_microbes(self, metabolite):
        ontology = metabolite.find("hmdb:ontology", self.namespace)
        for root in ontology.findall("hmdb:root", self.namespace):
            if self.get_text(root, "term") != "Disposition":
                continue
            for d in root.findall(".//hmdb:descendant", self.namespace):
                if self.get_text(d, "term") == "Microbe":
                    return sorted(
                        {
                            t.text.lower().strip()
                            for t in d.findall(".//hmdb:term", self.namespace)
                            if t.text and t.text.strip().lower() != "microbe"
                        }
                    )

    def get_diseases(self, metabolite):
        disease_names = set()
        diseases_elem = metabolite.find("hmdb:diseases", self.namespace)
        if diseases_elem is not None:
            for disease_elem in diseases_elem.findall("hmdb:disease", self.namespace):
                name_elem = disease_elem.find("hmdb:name", self.namespace)
                if name_elem is not None and name_elem.text:
                    disease_name = name_elem.text.strip().lower()
                    disease_names.add(disease_name)

        return sorted(disease_names)

    def get_experimental_properties(self, metabolite, prop_name):
        props = metabolite.find("hmdb:experimental_properties", self.namespace)
        if props is not None:
            for prop in props.findall("hmdb:property", self.namespace):
                kind = self.get_text(prop, "kind")
                if kind and kind.lower() == prop_name.lower():
                    return self.get_text(prop, "value")
        return None

    def get_molecular_weights(self, metabolite):
        avg = self.get_text(metabolite, "average_molecular_weight")
        mono = self.get_text(metabolite, "monisotopic_molecular_weight")

        weights = {}
        if avg:
            weights["average_molecular_weight"] = float(avg)
        if mono:
            weights["monoisotopic_molecular_weight"] = float(mono)
        return weights

    def get_primary_id(self, metabolite):
        # (tag, prefix) pairs for ID hierarchy
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

        def classify_kegg(val):
            if val.startswith("C"):
                return "KEGG.COMPOUND"
            elif val.startswith("G"):
                return "KEGG.GLYCAN"
            elif val.startswith("D"):
                return "KEGG.DRUG"
            return "KEGG"

        xrefs, primary_id = {}, None
        for tag, prefix in id_hierarchy:
            val = self.get_text(metabolite, tag)
            if not val:
                continue
            if tag == "kegg_id":
                prefix = classify_kegg(val)
            elif tag == "smiles":
                prefix = ""

            curie = f"{prefix}:{val}"
            if prefix == "PUBCHEM.COMPOUND":
                key = "pubchem_cid"
            elif prefix == "HMDB":
                key = "hmdb"
            elif prefix == "METLIN":
                key = "metlin"
            else:
                key = prefix.lower().split(".")[0]

            if primary_id is None:
                primary_id = curie

            xrefs.setdefault(key, curie)
        return primary_id, xrefs

    def get_references(self, description, microbes):
        reference_indicators = ["PMID", "DOI", "Wikipedia", "http", "www", ":"]
        matches = self.parenthetical_pattern.findall(description or "")
        microbes = microbes or []

        for sentence, paren in matches:
            if any(ind.lower() in paren.lower() for ind in reference_indicators):
                matched_microbes = [m for m in microbes if m.lower() in sentence.lower()]
                if not matched_microbes:
                    continue

                paren = re.sub(r"\bPMID[\s:]+(\d+)", r"PMID:\1", paren, flags=re.IGNORECASE)
                refs = [r.strip() for r in re.split(r"[;|,]", paren) if ":" in r and "CAS" not in r]

                if refs:
                    ref_dict = {}
                    for ref in refs:
                        try:
                            prefix, value = ref.split(":", 1)
                        except ValueError:
                            continue

                        key = (
                            "pmid"
                            if prefix.lower() == "pmid"
                            else "doi"
                            if prefix.lower() == "doi"
                            else "wikidata"
                            if "wikipedia" in prefix.lower()
                            else "url"
                            if "http" in prefix.lower() or "www" in prefix.lower()
                            else "article"
                        )

                        value = (
                            int(value.strip()) if key == "pmid" and value.strip().isdigit() else ref
                        )
                        ref_dict.setdefault(key, []).append(value)

                    for k in list(ref_dict):
                        if isinstance(ref_dict[k], list) and len(ref_dict[k]) == 1:
                            ref_dict[k] = ref_dict[k][0]

                    if "pmid" in ref_dict:
                        first_pmid = (
                            ref_dict["pmid"][0]
                            if isinstance(ref_dict["pmid"], list)
                            else ref_dict["pmid"]
                        )
                        ref_dict["id"] = f"PMID:{first_pmid}"
                    elif "url" in ref_dict:
                        first_url = (
                            ref_dict["url"][0]
                            if isinstance(ref_dict["url"], list)
                            else ref_dict["url"]
                        )
                        ref_dict["id"] = f"JournalArticle:{first_url}"
                    elif "doi" in ref_dict:
                        first_doi = (
                            ref_dict["doi"][0].split(":")[1].strip()
                            if isinstance(ref_dict["doi"], list)
                            else ref_dict["doi"].split(":")[1].strip()
                        )
                        ref_dict["id"] = f"doi:{first_doi}"
                    elif "wikidata" in ref_dict:
                        ref_dict["id"] = "Wikipedia"
                    elif "article" in ref_dict:
                        first_article = (
                            ref_dict["article"][0]
                            if isinstance(ref_dict["article"], list)
                            else ref_dict["article"]
                        )
                        ref_dict["id"] = first_article

                    ref_dict["type"] = "biolink:Publication"
                    return ref_dict

        return {}

    def get_anatomical_entities(self, metabolite):
        output = {
            "cellular_component": [],
            "biosample": [],
            "anatomical_entity": [],
        }

        bio_prop = metabolite.find("hmdb:biological_properties", self.namespace)
        if bio_prop is not None:
            cell_elem = bio_prop.find("hmdb:cellular_locations", self.namespace)
            if cell_elem is not None:
                for cell in cell_elem.findall("hmdb:cellular", self.namespace):
                    txt = cell.text.strip().lower()
                    if txt and txt.lower() not in output["cellular_component"]:
                        output["cellular_component"].append(txt.lower())

            specimen_elem = bio_prop.find("hmdb:biospecimen_locations", self.namespace)
            if specimen_elem is not None:
                for specimen in specimen_elem.findall("hmdb:biospecimen", self.namespace):
                    txt = specimen.text.strip()
                    if txt and txt.lower() not in output["biosample"]:
                        output["biosample"].append(txt.lower())

            tissue_elem = bio_prop.find("hmdb:tissue_locations", self.namespace)
            if tissue_elem is not None:
                for tissue in tissue_elem.findall("hmdb:tissue", self.namespace):
                    txt = tissue.text.strip()
                    if txt and txt.lower() not in output["anatomical_entity"]:
                        output["anatomical_entity"].append(txt.lower())

        return output

    def get_medi_references(self, disease_elem):
        pmids = []
        references_elem = disease_elem.find("hmdb:references", self.namespace)
        if references_elem is not None:
            for ref in references_elem.findall("hmdb:reference", self.namespace):
                pmid_elem = ref.find("hmdb:pubmed_id", self.namespace)
                if pmid_elem is not None and pmid_elem.text:
                    try:
                        pmids.append(int(pmid_elem.text.strip()))
                    except ValueError:
                        continue
        if pmids:
            pmids = sorted(set(pmids))
            pmid_value = pmids[0] if len(pmids) == 1 else pmids

            return {"id": f"PMID:{pmids[0]}", "pmid": pmid_value, "type": "biolink:Publication"}

    def build_metabolite_node(
        self,
        metabolite,
        primary_id: str,
        xrefs: dict,
        anatomical_entities: dict | None = None,
    ):
        name = self.get_text(metabolite, "name")
        state = self.get_text(metabolite, "state")

        synonyms_elem = metabolite.find("hmdb:synonyms", self.namespace)
        synonyms = self.get_list(synonyms_elem, "synonym") if synonyms_elem is not None else []

        node = {
            "id": primary_id,
            "name": name.lower() if name else None,
            "synonym": synonyms,
            "description": self.get_text(metabolite, "description"),
            "chemical_formula": self.get_text(metabolite, "chemical_formula"),
            "molecular_weight": self.get_molecular_weights(metabolite),
            "state": state.lower() if state else None,
            "water_solubility": self.get_experimental_properties(metabolite, "water_solubility"),
            "logp": self.get_experimental_properties(metabolite, "logp"),
            "melting_point": self.get_experimental_properties(metabolite, "melting_point"),
            "type": "biolink:SmallMolecule",
            "xrefs": xrefs,
        }

        if anatomical_entities:
            node.update(
                {
                    "cellular_component": anatomical_entities.get("cellular_component"),
                    "biosample": anatomical_entities.get("biosample"),
                    "anatomical_entity": anatomical_entities.get("anatomical_entity"),
                }
            )

        return self.remove_empty_none_values(node)

    def parse_microbe_metabolite(self):
        """Parse the HMDB XML for microbe-metabolite associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()
        if not self.cached_taxon_info:
            cache_metabolite_data(self.input_xml)

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            microbes = self.get_microbes(metabolite)
            description = self.get_text(metabolite, "description")
            references = self.get_references(description, microbes)

            association_node = {
                "predicate": "biolink:OrganismTaxonToChemicalEntityAssociation",
                "type": "has_metabolic_interaction_with",
                "primary_knowledge_source": "infores:hmdb_v5",
                "publication": references,
            }
            if "publication" in association_node:
                association_node["evidence_type"] = "ECO:0000305"  # manual assertion
            else:
                association_node["evidence_type"] = "ECO:0000000"  # unknown evidence
            association_node = self.remove_empty_none_values(association_node)

            primary_id, xrefs = self.get_primary_id(metabolite)
            object_node = self.build_metabolite_node(
                metabolite,
                primary_id,
                xrefs,
                anatomical_entities=self.get_anatomical_entities(metabolite),
            )

            if not microbes:
                continue
            for microbe in microbes:
                if microbe in self.cached_taxon_info:
                    subject_node = self.cached_taxon_info[microbe].copy()
                    subject_node["original_name"] = microbe.lower().strip()
                    subject_node["type"] = "biolink:OrganismTaxon"
                    subject_node["organism_type"] = get_organism_type(subject_node)
                    subject_node = self.remove_empty_none_values(subject_node)

                    yield {
                        "_id": str(uuid.uuid4()),
                        "association": association_node,
                        "object": object_node,
                        "subject": subject_node,
                    }

    def parse_metabolite_disease(self):
        """Parse the HMDB XML for metabolite-disease associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            subject_node = self.build_metabolite_node(
                metabolite,
                primary_id,
                xrefs,
                anatomical_entities=self.get_anatomical_entities(metabolite),
            )

            diseases = self.get_diseases(metabolite)
            diseases_elem = metabolite.find("hmdb:diseases", self.namespace)
            if diseases_elem is not None:
                for disease_elem in diseases_elem.findall("hmdb:disease", self.namespace):
                    references = self.get_medi_references(disease_elem)
                    if references is None:
                        continue

                    association_node = {
                        "id": "RO:0000087",  # corresponds to has role
                        "predicate": "biolink:ChemicalToDiseaseOrPhenotypicFeatureAssociation",
                        "type": "has_role_in",
                        "primary_knowledge_source": "infores:hmdb_v5",
                        "evidence_type": "ECO:0000305",  # manual assertion
                        "publication": references,
                    }
                    association_node = self.remove_empty_none_values(association_node)

                    for disease in diseases:
                        disease_key = disease.strip().lower()
                        if disease_key in self.cached_disease_info:
                            object_node = self.cached_disease_info[disease_key].copy()
                            object_node = self.remove_empty_none_values(object_node)

                            yield {
                                "_id": str(uuid.uuid4()),
                                "association": association_node,
                                "object": object_node,
                                "subject": subject_node,
                            }

    def parse_metabolite_protein(self):
        """Parse the HMDB XML for metabolite-protein associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            subject_node = self.build_metabolite_node(
                metabolite,
                primary_id,
                xrefs,
                anatomical_entities=self.get_anatomical_entities(metabolite),
            )

            prot_elem = metabolite.find("hmdb:protein_associations", self.namespace)
            if prot_elem is None:
                continue
            for protein in prot_elem.findall("hmdb:protein", self.namespace):
                object_node = {
                    "id": None,
                    "name": self.get_text(protein, "gene_name"),
                    "full_name": self.get_text(protein, "name"),
                    "description": None,
                    "type": "biolink:Protein",
                    "protein_type": (ptype := self.get_text(protein, "protein_type"))
                    and ptype.lower(),
                    "xrefs": {},
                }

                prot_accession = self.get_text(protein, "protein_accession")
                object_node["xrefs"]["hmdbp"] = f"HMDBP:{prot_accession}"

                uniprot_id = protein.findtext(
                    "hmdb:uniprot_id", default="", namespaces=self.namespace
                ).strip()
                if uniprot_id:
                    object_node["id"] = f"UniProtKB:{uniprot_id}"
                    object_node["xrefs"]["uniprotkb"] = f"UniProtKB:{uniprot_id}"
                    for d in self.cached_protein_function:
                        if uniprot_id in d:
                            object_node["description"] = d[uniprot_id]["description"]
                elif prot_accession:
                    object_node["id"] = f"HMDBP:{prot_accession}"

                object_node = self.remove_empty_none_values(object_node)

                association_node = {
                    "id": "RO:0002434",
                    "predicate": "biolink:ChemicalGeneInteractionAssociation",
                    "type": "interacts_with",
                    "primary_knowledge_source": "infores:hmdb_v5",
                    "evidence_type": "ECO:0000000",  # unknown evidence
                }

                yield {
                    "_id": str(uuid.uuid4()),
                    "association": association_node,
                    "object": object_node,
                    "subject": subject_node,
                }

    def parse_metabolite_pathway(self):
        """Parse the HMDB XML for metabolite-pathway associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            subject_node = self.build_metabolite_node(
                metabolite,
                primary_id,
                xrefs,
                anatomical_entities=self.get_anatomical_entities(metabolite),
            )

            smpdb_pw = self.cached_pathway_descr
            bio_prop = metabolite.find("hmdb:biological_properties", self.namespace)
            pathways_elem = bio_prop.find("hmdb:pathways", self.namespace)
            if pathways_elem is None:
                continue
            for pw in pathways_elem.findall("hmdb:pathway", self.namespace):
                pw_name = self.get_text(pw, "name")
                kegg_map = self.get_text(pw, "kegg_map_id")
                cache = smpdb_pw.get(pw_name.lower())
                smpdb_id = cache["smpdb"] if cache else None
                descr = cache["description"] if cache else None

                object_node = {
                    "id": smpdb_id or (f"KEGG:{kegg_map}" if kegg_map else None),
                    "name": pw_name.lower(),
                    "description": descr,
                    "type": "biolink:Pathway",
                    "xrefs": {
                        "smpdb": smpdb_id,
                        "kegg": f"KEGG:{kegg_map}" if kegg_map else None,
                    },
                }
                object_node = self.remove_empty_none_values(object_node)

                association_node = {
                    "id": "RO:0000056",
                    "predicate": "biolink:ChemicalToPathwayAssociation",
                    "type": "participates_in",
                    "primary_knowledge_source": "infores:hmdb_v5",
                    "evidence_type": "ECO:0000000",  # unknown evidence
                }

                _id = (
                    f"{subject_node['id'].split(':')[1]}_participates_in_{object_node['id'].split(':')[1]}"
                    if "id" in object_node and "id" in subject_node
                    else str(uuid.uuid4())
                )

                yield {
                    "_id": _id,
                    "association": association_node,
                    "object": object_node,
                    "subject": subject_node,
                }


class HMDB_Protein_Parse(XMLParseHelper):
    def __init__(self, input_xml):
        super().__init__(input_xml)
        self.input_xml = input_xml
        self.uniprot2entrezgene = load_pickle("hmdbp_uniprot2entrezgene.pkl")
        self.protein_func = {
            uniprot: info
            for d in load_pickle("hmdbp_uniprot_protein_functions.pkl")
            for uniprot, info in d.items()
        }
        self.gene_summary = {
            str(entrezgene): info
            for d in load_pickle("hmdbp_entrezgene_summaries.pkl")
            for entrezgene, info in d.items()
        }
        self.cached_pathway_descr = load_pickle("smpdb_pathway_descriptions.pkl")
        self.cached_go_descr = load_pickle("hmdbp_go_definitions.pkl")

    def get_list_of_region_list(self, elem, tag):
        ranges = []
        for e in elem.findall(f"hmdb:{tag}", self.namespace):
            if not e.text or "-" not in e.text:
                continue
            start, end = (part.strip() for part in e.text.split("-", 1))
            if start.isdigit() and end.isdigit():
                ranges.append([int(start), int(end)])
                continue
        return ranges

    def get_protein_properties(self, protein):
        props = protein.find("hmdb:protein_properties", self.namespace)
        if props is not None:
            residue_num = self.get_text(props, "residue_number")
            molecular_weight = self.get_text(props, "molecular_weight")
            pi = self.get_text(props, "theoretical_pi")
            tm_elem = props.find("hmdb:transmembrane_regions", self.namespace)
            tm_regions = (
                self.get_list_of_region_list(tm_elem, "region") if tm_elem is not None else []
            )
            sig_elem = props.find("hmdb:signal_regions", self.namespace)
            sig_regions = (
                self.get_list_of_region_list(sig_elem, "region") if sig_elem is not None else []
            )
            prot_seq = self.get_text(props, "polypeptide_sequence")
            pfams_elem = props.find("hmdb:pfams", self.namespace)
            pfam_list = [
                {
                    "id": f"PFAM:{pfam_id}",
                    "name": name.lower() if name else None,
                }
                for pfam in pfams_elem.findall("hmdb:pfam", self.namespace)
                if (pfam_id := self.get_text(pfam, "pfam_id"))
                and (name := self.get_text(pfam, "name"))
            ]
            return {
                "residue_num": int(residue_num) if residue_num else None,
                "molecular_weight": float(molecular_weight) if molecular_weight else None,
                "theoretical_pi": float(pi) if pi else None,
                "transmembrane_region": tm_regions if tm_regions else None,
                "signal_region": sig_regions if sig_regions else None,
                "protein_seq": prot_seq if prot_seq else None,
                "pfam": pfam_list if pfam_list else None,
            }

    def get_gene_properties(self, protein):
        props = protein.find("hmdb:gene_properties", self.namespace)
        if props is not None:
            chrom_loc_raw = self.get_text(props, "chromosome_location")
            locus = self.get_text(props, "locus")
            gene_seq = self.get_text(props, "gene_sequence")

            chrom_loc = None
            if chrom_loc_raw:
                candidate = (
                    chrom_loc_raw.split(":", 1)[-1] if ":" in chrom_loc_raw else chrom_loc_raw
                )
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

    def get_primary_id(self, protein):
        """
        Extract primary ID with hiearchy. pfams and pdbs are not included in the hierarchy yet.
        :param protein:
        :return:
        """
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
            val = self.get_text(protein, tag)
            if not val:
                continue

            curie = f"{prefix}:{val}"
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

    def get_references(self, protein):
        pmids = []
        ref_elem = protein.find("hmdb:general_references", self.namespace)
        if ref_elem is not None:
            for ref in ref_elem.findall("hmdb:reference", self.namespace):
                pmid_elem = ref.find("hmdb:pubmed_id", self.namespace)
                if pmid_elem is not None and pmid_elem.text:
                    pmids.append(int(pmid_elem.text.strip()))

        if pmids:
            pmids = sorted(set(pmids))
            pmid_value = pmids[0] if len(pmids) == 1 else pmids
            return {"id": f"PMID:{pmids[0]}", "pmid": pmid_value, "type": "biolink:Publication"}

    def lookup_uniprot(self, uniprot_id: str) -> tuple[str | None, str | None]:
        descr = self.protein_func.get(uniprot_id, {}).get("description")
        entrez_map = self.uniprot2entrezgene.get(uniprot_id, {})
        entrez_curie = entrez_map.get("gene_id")
        return descr, entrez_curie

    def lookup_entrez(self, entrez_curie: str | None) -> tuple[str | None, str | None]:
        if not entrez_curie or ":" not in entrez_curie:
            return None, None
        entrez_id = entrez_curie.split(":", 1)[1]
        gene_descr = self.gene_summary.get(entrez_id, {}).get("description")
        return entrez_id, gene_descr

    def build_protein_node(self, protein):
        name = self.get_text(protein, "gene_name")
        full_name = self.get_text(protein, "name")
        synonyms_elem = protein.find("hmdb:synonyms", self.namespace)
        synonyms = self.get_list(synonyms_elem, "synonym") if synonyms_elem is not None else []
        primary_id, xrefs = self.get_primary_id(protein)
        prot_func = self.get_text(protein, "general_function")
        prot_spec_func = self.get_text(protein, "specific_function")
        prot_type = self.get_text(protein, "protein_type")
        cellular_com_elem = protein.find("hmdb:subcellular_locations", self.namespace)
        cellular_components = (
            self.get_list(cellular_com_elem, "subcellular_location")
            if cellular_com_elem is not None
            else []
        )
        pdb_elem = protein.find("hmdb:pdb_ids", self.namespace)
        pdbs = self.get_list(pdb_elem, "pdb_id") if pdb_elem is not None else []
        xrefs["pdb"] = pdbs if pdb_elem is not None else []

        prot_props = self.get_protein_properties(protein)
        xrefs["pfam"] = prot_props["pfam"] if "pfam" in prot_props else []
        gene_props = self.get_gene_properties(protein)

        uniprot_id = self.get_text(protein, "uniprot_id")
        if not uniprot_id:
            return None
        prot_descr, entrez_curie = self.lookup_uniprot(uniprot_id)
        entrez_id, gene_descr = self.lookup_entrez(entrez_curie)

        if entrez_id:
            xrefs["entrezgene"] = f"NCBIGene:{entrez_id}"

        protein_node = {
            "id": primary_id,
            "name": name if name else None,
            "full_name": full_name.lower() if full_name else None,
            "synonym": synonyms,
            "description": prot_descr,
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
            "gene_description": gene_descr,
            "type": "biolink:Protein",
            "protein_type": prot_type.lower() if prot_type else None,
            "cellular_component": cellular_components,
            "xrefs": xrefs,
        }
        return self.remove_empty_none_values(protein_node)

    def parse_protein_pathway(self):
        """Parse the HMDB XML for protein-pathway associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()
        if not self.uniprot2entrezgene:
            cache_protein_data(self.input_xml)

        for protein in root.findall("hmdb:protein", self.namespace):
            subject_node = self.build_protein_node(protein)

            smpdb_pw = self.cached_pathway_descr
            pathways_elem = protein.find("hmdb:pathways", self.namespace)
            if pathways_elem is None:
                continue
            for pw in pathways_elem.findall("hmdb:pathway", self.namespace):
                pw_name = self.get_text(pw, "name")
                kegg_map = self.get_text(pw, "kegg_map_id")
                cache = smpdb_pw.get(pw_name.lower())
                smpdb_id = cache["smpdb"] if cache else None
                descr = cache["description"] if cache else None

                object_node = {
                    "id": smpdb_id or (f"KEGG:{kegg_map}" if kegg_map else None),
                    "name": pw_name.lower(),
                    "description": descr,
                    "type": "biolink:Pathway",
                    "xrefs": {
                        "smpdb": smpdb_id,
                        "kegg": f"KEGG:{kegg_map}" if kegg_map else None,
                    },
                }
                object_node = self.remove_empty_none_values(object_node)

                publication = self.get_references(protein)
                association_node = {
                    "id": "RO:0000056",
                    "predicate": "biolink:GeneToPathwayAssociation",
                    "type": "participates_in",
                    "primary_knowledge_source": "infores:hmdb_v5",
                    "evidence_type": "ECO:0000305",  # manual assertion
                    "publication": publication,
                }
                association_node = self.remove_empty_none_values(association_node)

                _id = (
                    f"{subject_node['id'].split(':')[1]}_participates_in_{object_node['id'].split(':')[1]}"
                    if "id" in object_node and "id" in subject_node
                    else str(uuid.uuid4())
                )

                yield {
                    "_id": _id,
                    "association": association_node,
                    "object": object_node,
                    "subject": subject_node,
                }

    def parse_protein_biological_process(self):
        """Parse the HMDB XML for protein-biological process associations."""
        tree = ET.parse(self.input_xml)
        root = tree.getroot()

        for protein in root.findall("hmdb:protein", self.namespace):
            subject_node = self.build_protein_node(protein)
            publication = self.get_references(protein)

            bp_elem = protein.find("hmdb:go_classifications", self.namespace)
            if bp_elem is None:
                continue
            for bp in bp_elem.findall("hmdb:go_class", self.namespace):
                if self.get_text(bp, "category") != "Biological process":
                    continue

                go_name = self.get_text(bp, "description")
                go_id = self.get_text(bp, "go_id")
                go_descr = self.cached_go_descr.get(go_id, {}).get("description") if go_id else None

                object_node = {
                    "id": go_id,
                    "name": go_name.lower() if go_name else None,
                    "description": go_descr,
                    "type": "biolink:BiologicalProcess",
                    "xrefs": {"go": go_id} if go_id else {},
                }
                object_node = self.remove_empty_none_values(object_node)

                association_node = {
                    "id": "RO:0002331",
                    "predicate": "biolink:MacromolecularMachineToBiologicalProcessAssociation",
                    "type": "involved_in",
                    "primary_knowledge_source": "infores:hmdb_v5",
                    "evidence_type": "ECO:0000305",  # manual assertion
                    "publication": publication,
                }
                association_node = self.remove_empty_none_values(association_node)

                _id = (
                    f"{subject_node['id'].split(':', 1)[1]}_involved_in_{object_node['id'].split(':', 1)[1]}"
                    if "id" in object_node and "id" in subject_node
                    else str(uuid.uuid4())
                )
                yield {
                    "_id": _id,
                    "association": association_node,
                    "object": object_node,
                    "subject": subject_node,
                }


def load_hmdb_data(data_dir="downloads"):
    """
    Load HMDB all association data including:
    - Microbe-Metabolite
    - Metabolite-Disease, Metabolite-Protein, Metabolite-Pathway
    - Protein-Pathway, Protein-Biological Process
    """
    data_dir = os.path.abspath(data_dir)
    metabolite_xml = os.path.join(data_dir, "hmdb_metabolites.xml")
    if not os.path.isfile(metabolite_xml):
        metabolite_zip = os.path.join(data_dir, "hmdb_metabolites.zip")
        metabolite_xml = extract_file_from_zip(metabolite_zip, "hmdb_metabolites.xml")

    protein_xml = os.path.join(data_dir, "hmdb_proteins.xml")
    if not os.path.isfile(protein_xml):
        protein_zip = os.path.join(data_dir, "hmdb_proteins.zip")
        protein_xml = extract_file_from_zip(protein_zip, "hmdb_proteins.xml")

    metabolite_parser = HMDB_Metabolite_Parse(metabolite_xml)
    protein_parser = HMDB_Protein_Parse(protein_xml)

    yield from chain(
        metabolite_parser.parse_microbe_metabolite(),
        metabolite_parser.parse_metabolite_disease(),
        metabolite_parser.parse_metabolite_protein(),
        metabolite_parser.parse_metabolite_pathway(),
        protein_parser.parse_protein_pathway(),
        protein_parser.parse_protein_biological_process(),
    )


def cache_hmdb_db(data_dir="downloads"):
    """Cache HMDB data for faster access in the future."""
    print(f"--- Caching HMDB data in {data_dir}... ---")
    data_dir = os.path.abspath(data_dir)
    metabolite_xml = os.path.join(data_dir, "hmdb_metabolites.xml")
    if not os.path.isfile(metabolite_xml):
        metabolite_zip = os.path.join(data_dir, "hmdb_metabolites.zip")
        metabolite_xml = extract_file_from_zip(metabolite_zip, "hmdb_metabolites.xml")

    protein_xml = os.path.join(data_dir, "hmdb_proteins.xml")
    if not os.path.isfile(protein_xml):
        protein_zip = os.path.join(data_dir, "hmdb_proteins.zip")
        protein_xml = extract_file_from_zip(protein_zip, "hmdb_proteins.xml")

    metabolite_parser = HMDB_Metabolite_Parse(metabolite_xml)
    protein_parser = HMDB_Protein_Parse(protein_xml)

    hmdb_combined = {
        "microbe-metabolite": list(tqdm(metabolite_parser.parse_microbe_metabolite())),
        "metabolite-disease": list(tqdm(metabolite_parser.parse_metabolite_disease())),
        "metabolite-protein": list(tqdm(metabolite_parser.parse_metabolite_protein())),
        "metabolite-pathway": list(tqdm(metabolite_parser.parse_metabolite_pathway())),
        "protein-biological_process": list(tqdm(protein_parser.parse_protein_biological_process())),
        "protein-pathway": list(tqdm(protein_parser.parse_protein_pathway())),
    }

    print("--- Saving parsed records to pickle file... ---")
    save_pickle(
        hmdb_combined,
        "hmdb_v5_parsed_records.pkl",
    )
    print("--- Caching complete. ---")


if __name__ == "__main__":
    start = time.time()
    cache_hmdb_db()
    # cache_metabolite_data(os.path.join("downloads", "hmdb_metabolites.xml"))
    # cache_protein_data(os.path.join("downloads", "hmdb_proteins.xml"))
    # recs = [rec for rec in load_hmdb_data()]
    # for rec in recs:
    #     print(rec)
    end = time.time()
    print(f"Total time: {(end - start)/60:.2f} minutes.")  # it took ~28 minutes with cache
