import os
import pathlib
import pickle
import ssl
import time
import zipfile
from typing import Iterator

import biothings_client as bt
import requests
import text2term
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from lxml import etree as ET

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)
load_dotenv()


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


def extract_xml_from_zip(zip_path, expected_filename="hmdb_metabolites.xml"):
    extract_dir = os.path.dirname(zip_path)
    extracted_path = os.path.join(extract_dir, expected_filename)

    if not os.path.exists(extracted_path):
        with zipfile.ZipFile(zip_path, "r") as zip_f:
            zip_f.extract(expected_filename, path=extract_dir)

    return extracted_path


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
    taxon_names: list[str],
    email,
    sleep=0.34,
) -> dict:
    """Map taxonomy names to NCBI taxonomy ids using entrez API
    Entrez is good at mapping recent reclassified taxonomy names that are outdated
    but very slow due to no batch query allowed, so expensive to query, recommend cache the output

    :param taxon_names:
    :param email:
    :param sleep:
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
          "id":"taxid:28450",
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
                "id": f"taxid:{int(info['_id'])}",
                "taxid": int(info["_id"]),
                "name": info["scientific_name"],
                "parent_taxid": int(info["parent_taxid"]),
                "lineage": info["lineage"],
                "rank": info["rank"],
            }
        else:
            taxon[info["query"]] = {
                "id": f"taxid:{int(info['query'])}",
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
    API_KEY = os.getenv("NCIT_API_KEY")
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
    '28450': {'id': 'taxid:28450', 'taxid': 28450, 'name': 'burkholderia pseudomallei', 'parent_taxid': 111527, 'lineage': [28450, 111527, 32008, 119060, 80840, 28216, 1224, 3379134, 2, 131567, 1], 'rank': 'species', 'description': 'A species of aerobic, Gram-negative, rod shaped bacteria assigned to the phylum Proteobacteria.
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


def cache_data(input_xml):
    microbe_names = get_all_microbe_names(input_xml)
    save_pickle(list(set(microbe_names)), "hmdb_v5_microbe_names.pkl")

    microbes4query = load_pickle("hmdb_v5_microbe_names.pkl")
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


class HMDBParse:
    def __init__(self, input_xml):
        self.namespace = {"hmdb": "http://www.hmdb.ca"}
        self.input_xml = input_xml
        self.cached_taxon_info = load_pickle("original_taxon_name2taxid.pkl")

    def get_text(self, elem, tag):
        child = elem.find(f"hmdb:{tag}", self.namespace)
        return child.text.strip() if child and child.text else None

    def get_list(self, elem, tag):
        return [e.text.lower() for e in elem.findall(f"hmdb:{tag}", self.namespace) if e.text]

    def get_experimental_properties(self, metabolite, prop_name):
        props = metabolite.find("hmdb:experimentalProperties", self.namespace)
        if props:
            for prop in props.findall("hmdb:property", self.namespace):
                kind = self.get_text(prop, "kind")
                if kind and kind.lower() == prop_name.lower():
                    return self.get_text(prop, "value")
        return None

    def get_molecular_weights(self, metabolite):
        return {
            "average_molecular_weight": self.get_text(metabolite, "average_molecular_weight"),
            "monisotopic_molecular_weight": self.get_text(
                metabolite, "monisotopic_molecular_weight"
            ),
        }

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

    def get_primary_id(self, metabolite):
        id_hierarchy = [
            ("pubchem_compound_id", "PUBCHEM.COMPOUND"),
            ("inchikey", "INCHIKEY"),
            ("drugbank_id", "DRUGBANK"),
            ("chebi_id", "CHEBI"),
            ("chembl_id", "CHEMBL.COMPOUND"),
            ("accession", "HMDB"),
            ("cas_registry_number", "CAS"),
            ("kegg_id", None),
            ("chemspider_id", "chemspider"),
            ("foodb_id", "foodb.compound"),
            ("bigg_id", "BIGG.METABOLITE"),
        ]

        def classify_kegg(val):
            if val.startswith("C"):
                return "KEGG.COMPOUND"
            elif val.startswith("G"):
                return "KEGG.GLYCAN"
            elif val.startswith("D"):
                return "KEGG.DRUG"
            return "KEGG"

        xrefs = {}
        primary_id = None

        for tag, prefix in id_hierarchy:
            val = self.get_text(metabolite, tag)
            if val:
                if tag == "kegg_id":
                    prefix = classify_kegg(val)
                curie = f"{prefix}:{val}" if prefix != "HMDB" else val.replace("HMDB", "")
                if not primary_id:
                    primary_id = curie
                else:
                    xrefs[prefix] = val
        return primary_id, xrefs

    def parse(self):
        tree = ET.parse(self.input_xml)
        root = tree.getroot()
        if not os.path.exists("original_taxon_name2taxid.pkl"):
            cache_data(self.input_xml)

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            rec = {
                "_id": None,
                "association": {},
                "object": {},
                "subject": {},
            }

            object_node = {
                "id": primary_id,
                "name": self.get_text(metabolite, "name").lower()
                if self.get_text(metabolite, "name")
                else None,
                "synonym": self.get_list(
                    metabolite.find("hmdb:synonyms", self.namespace), "synonym"
                )
                if metabolite.find("hmdb:synonyms", self.namespace)
                else [],
                "description": self.get_text(metabolite, "description"),
                "chemical_formula": self.get_text(metabolite, "chemical_formula"),
                "molecular_weight": self.get_molecular_weights(metabolite),
                "state": self.get_experimental_properties(metabolite, "state"),
                "water_solubility": self.get_experimental_properties(
                    metabolite, "water_solubility"
                ),
                "logp": self.get_experimental_properties(metabolite, "logp"),
                "melting_point": self.get_experimental_properties(metabolite, "melting_point"),
                "type": "biolink:SmallMolecule",
                "xrefs": xrefs,
            }
            rec["object"] = object_node

            microbes = self.get_microbes(metabolite)
            taxon_info = self.cached_taxon_info
            if not microbes:
                continue
            for microbe in microbes:
                if microbe in taxon_info:
                    subject_node = taxon_info[microbe]
                    subject_node["original_name"] = microbe.lower().strip()
                    subject_node["type"] = "biolink:OrganismTaxon"
                    rec["subject"] = subject_node

                    _id = f"{rec['subject']['id'].split(':')[1]}_associated_with_{rec['object']['id'].split(':')[1]}"
                    rec["_id"] = _id
                    yield rec


if __name__ == "__main__":
    zip_path = os.path.join("downloads", "hmdb_metabolites.zip")
    hmdb_xml = extract_xml_from_zip(zip_path)
    parser = HMDBParse(xml_path)
    for rec in parser.parse():
        print(rec)
