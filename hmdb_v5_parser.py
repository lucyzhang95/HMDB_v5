import os
import ssl
import zipfile

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
    idx = tag.rfind("}")
    # rfind() method "not found" == -1
    if idx != -1:  # if idx is not "not found"
        tag = tag[idx + 1 :]
    return tag


def ete3_taxon_name2taxid(taxon_names: list) -> dict:
    """Use ete3 to map taxonomy names to NCBI taxonomy ids
    ete3 is good at mapping exact taxonomy names and fast without accessing API

    :param taxon_names: a list of taxon names (the values of the output of preprocess_taxon_name)
    :return: a dictionary mapping taxon names to taxid numbers
    {'human papillomavirus 11': 10580, 'veillonella sp.': 1926307, ...}
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
            ete3_mapped[name] = {"taxid": int(taxid[0])}
    return ete3_mapped
