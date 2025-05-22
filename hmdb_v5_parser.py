import os
import pathlib
import pickle
import ssl
import time
import zipfile
from typing import Iterator

import requests
import text2term
from Bio import Entrez
from ete3 import NCBITaxa
from lxml import etree as ET

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    """
    :param obj:
    :param f_name: files should only be existing in the cache directory
    :return:
    """
    f_path = os.path.join(CACHE_DIR, f_name)
    with open(f_path, "wb") as in_f:
        pickle.dump(obj, in_f)


def load_pickle(f_name):
    f_path = os.path.join(CACHE_DIR, f_name)
    if os.path.exists(f_path):
        with open(f_path, "rb") as in_f:
            return pickle.load(in_f)
    return None


def strip_tag_namespace(tag: str) -> str:
    """Strip namespace from main element tag
    Each element tag from HMDB contains "http://www.hmdb.ca" as namespace
    e.g. {http://www.hmdb.ca}accession
    Remove the namespace for easier data processing

    :param tag: element tags (e.g. <accession>HMDB0000001</accession>)
    :return: original tag without namespace
    """
    idx = tag.rfind("}")  # rfind() "not found" == -1
    if idx != -1:  # if idx is not "not found"
        tag = tag[idx + 1 :]
    return tag


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

    for _, elem in ET.iterparse(str(input_xml), events=("end",), tag=f"{{{ns_uri}}}metabolite"):
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
    taxon_names: list[str], sleep=0.34, email="bazhang@scripps.edu"
) -> dict:
    """Map taxonomy names to NCBI taxonomy ids using entrez API
    Entrez is good at mapping recent reclassified taxonomy names that are outdated
    but very slow due to no batch query allowed, so expensive to query, recommend cache the output

    :param taxon_names:
    :param sleep:
    :param email:
    :return: a dictionary mapping taxon names to taxids,
    e.g., {'cellulomonas galba': {'taxid': 401861, 'mapping_tool': 'entrez'}, ...}
    """
    Entrez.email = email
    entrez_mapped = {}

    for name in set(taxon_names):
        try:
            handle = Entrez.esearch(db="taxonomy", term=name, retmode="xml", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                taxid = int(record["IdList"][0])
                entrez_mapped[name] = {"taxid": taxid, "mapping_tool": "entrez"}
        except Exception as e:
            print(f"Entrez query failed for '{name}': {e}")
        time.sleep(sleep)
    return entrez_mapped


def text2term_taxon_name2taxid(taxon_names: list[str], min_score=0.8) -> dict:
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


def manual_taxon_name2taxid():
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
        name: {"taxid": taxid, "mapping_tool": "manual"} for name, taxid in raw_mapping.items()
    }
    return manual_mapped


def get_ncit_taxon_description(taxon_names):
    API_KEY = "efd61c1d-74a2-4877-b4ff-37ba827a96bc"
    search_url = "https://data.bioontology.org/search"
    taxon_names = set(taxon_names)
    mapping_result = {}
    for name in taxon_names:
        params = {
            "q": name,
            "ontologies": "NCIT",
            "apikey": API_KEY,
        }
        response = requests.get(search_url, params=params)
        data = response.json()
        for result in data.get("collection", []):
            if result:
                ncit_output = {
                    "ncit": result.get("@id").split("#")[1],
                    "name": result.get("prefLabel").lower(),
                    "description": f"{result.get('definition')[0]} [NCIT]"
                    if "definition" in result
                    else "",
                }
                if ncit_output["name"] == name:
                    mapping_result[name] = ncit_output
                    del mapping_result[name]["name"]
    return mapping_result


if __name__ == "__main__":
    input_xml = os.path.join("downloads", "hmdb_metabolites.xml")
    microbe_names = get_all_microbe_names(input_xml)
    # save_pickle(list(set(microbe_names)), "hmdb_v5_microbe_names.pkl")
    microbes4query = [obj for obj in microbe_names]
    # ete3_mapped = ete3_taxon_name2taxid(microbes4query)
    # save_pickle(ete3_mapped, "ete3_name2taxid.pkl")
    ete3_cached = load_pickle("ete3_name2taxid.pkl")

    # no_hits = [name for name in set(microbe_names) if name not in ete3_mapped]
    microbe_cached = load_pickle("hmdb_v5_microbe_names.pkl")
    no_hits = [name for name in microbe_cached if name not in ete3_cached]
    # entrez_mapped = entrez_taxon_name2taxid(no_hits)
    # save_pickle(entrez_mapped, "entrez_name2taxid.pkl")
    entrez_cached = load_pickle("entrez_name2taxid.pkl")

    # no_hits2 = [
    #     name for name in set(microbe_names) if name not in ete3_mapped and name not in entrez_mapped
    # ]
    no_hits2 = [
        name for name in microbe_cached if name not in ete3_cached and name not in entrez_cached
    ]

    # text2term_mapped = text2term_taxon_name2taxid(no_hits2)
    # text2term_mapped = manual_correct_text2term_map(text2term_mapped)
    # save_pickle(text2term_mapped, "text2term_name2taxid.pkl")
    text2term_cached = load_pickle("text2term_name2taxid.pkl")

    # no_hits3 = [name for name in set(microbe_names) if name not in ete3_mapped and name not in entrez_mapped and name not in text2term_mapped]
    no_hits3 = [
        name
        for name in microbe_cached
        if name not in ete3_cached and name not in entrez_cached and name not in text2term_cached
    ]
    print(no_hits3)
