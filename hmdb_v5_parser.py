import os
import pathlib
import ssl
import time
import zipfile
from typing import Iterator

import requests
from Bio import Entrez
from ete3 import NCBITaxa
from lxml import etree as ET

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


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
