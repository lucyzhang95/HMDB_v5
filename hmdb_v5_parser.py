import asyncio
import os
import pathlib
import pickle
import re
import ssl
import time
import uuid
import zipfile
from typing import Dict, Iterator, List, Optional, Union

import aiohttp
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


def cache_data(input_xml):
    microbe_names = get_all_microbe_names(input_xml)
    save_pickle(list(set(microbe_names)), "hmdb_v5_microbe_names.pkl")

    # cache mapped taxon
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

    # cache mapped diseases
    diseases = get_all_diseases(hmdb_xml)
    save_pickle(diseases, "hmdb_v5_diseases.pkl")

    diseases = load_pickle("hmdb_v5_diseases.pkl")
    disease4text2term = [di_name for di_name, omim in diseases.items() if omim is None]
    text2term_mapped = text2term_disease_name2id(disease4text2term)
    save_pickle(text2term_mapped, "text2term_disease_name2id.pkl")
    disease2cui = [di_name for di_name in disease4text2term if di_name not in text2term_mapped]
    umls_mapped = get_cuis_sync(os.getenv("UMLS_API_KEY"), disease2cui)
    save_pickle(umls_mapped, "umls_disease_name2id.pkl")
    di_no_hit = [name for name, mapped in umls_mapped.items() if mapped["id"] == ""]
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

    di2id_path = os.path.join("cache", "all_disease_name2id.pkl")
    if os.path.exists(di2id_path):
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


class HMDBParse:
    def __init__(self, input_xml):
        self.namespace = {"hmdb": "http://www.hmdb.ca"}
        self.input_xml = input_xml
        self.cached_taxon_info = load_pickle("original_taxon_name2taxid.pkl")
        self.cached_disease_info = load_pickle("original_disease_name2id.pkl")
        self.parenthetical_pattern = re.compile(r"([^.?!]*?)\s*\(([^)]*?)\)")

    def get_text(self, elem, tag):
        child = elem.find(f"hmdb:{tag}", self.namespace)
        return child.text.strip() if child is not None and child.text else None

    def get_list(self, elem, tag):
        return [e.text.lower() for e in elem.findall(f"hmdb:{tag}", self.namespace) if e.text]

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

        xrefs = {}
        primary_id = None
        for tag, prefix in id_hierarchy:
            val = self.get_text(metabolite, tag)
            if val:
                if tag == "kegg_id":
                    prefix = classify_kegg(val)

                curie = f"{prefix}:{val}"
                if not primary_id:
                    primary_id = curie
                else:
                    key = (
                        "foodb"
                        if prefix == "foodb.compound"
                        else "hmdb"
                        if prefix == "HMDB"
                        else prefix.lower().split(".")[0]
                    )
                    xrefs[key] = curie

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

    def remove_empty_none_values(self, d):
        if isinstance(d, dict):
            return {
                k: self.remove_empty_none_values(v) for k, v in d.items() if v not in (None, {}, [])
            }
        elif isinstance(d, list):
            return [self.remove_empty_none_values(v) for v in d if v not in (None, {}, [])]
        else:
            return d

    def parse_mime(self):
        tree = ET.parse(self.input_xml)
        root = tree.getroot()
        cached_taxon_info = os.path.join("cache", "original_taxon_name2taxid.pkl")
        if not os.path.exists(cached_taxon_info):
            cache_data(self.input_xml)

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            microbes = self.get_microbes(metabolite)
            description = self.get_text(metabolite, "description")
            references = self.get_references(description, microbes)

            association_node = {
                "predicate": "biolink:OrganismTaxonToChemicalEntityAssociation",
                "infores": "hmdb_v5",
                "publication": references,
            }
            association_node = self.remove_empty_none_values(association_node)

            name = self.get_text(metabolite, "name")
            synonyms_elem = metabolite.find("hmdb:synonyms", self.namespace)
            logp = self.get_experimental_properties(metabolite, "logp")
            state = self.get_text(metabolite, "state")
            object_node = {
                "id": primary_id,
                "name": name.lower() if name else None,
                "synonym": self.get_list(synonyms_elem, "synonym")
                if synonyms_elem is not None
                else []
                if metabolite.find("hmdb:synonyms", self.namespace)
                else [],
                "description": description,
                "chemical_formula": self.get_text(metabolite, "chemical_formula"),
                "molecular_weight": self.get_molecular_weights(metabolite),
                "state": state.lower() if state else None,
                "water_solubility": self.get_experimental_properties(
                    metabolite, "water_solubility"
                ),
                "logp": logp,
                "melting_point": self.get_experimental_properties(metabolite, "melting_point"),
                "type": "biolink:SmallMolecule",
                "xrefs": xrefs,
            }
            object_node = self.remove_empty_none_values(object_node)

            if not microbes:
                continue
            for microbe in microbes:
                if microbe in self.cached_taxon_info:
                    subject_node = self.cached_taxon_info[microbe].copy()
                    subject_node["original_name"] = microbe.lower().strip()
                    subject_node["type"] = "biolink:OrganismTaxon"
                    subject_node = self.remove_empty_none_values(subject_node)

                    rec = {
                        "_id": str(uuid.uuid4()),
                        "association": association_node,
                        "object": object_node,
                        "subject": subject_node,
                    }

                    yield rec

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

    def parse_medi(self):
        tree = ET.parse(self.input_xml)
        root = tree.getroot()

        for metabolite in root.findall("hmdb:metabolite", self.namespace):
            primary_id, xrefs = self.get_primary_id(metabolite)
            description = self.get_text(metabolite, "description")
            diseases = self.get_diseases(metabolite)

            name = self.get_text(metabolite, "name")
            synonyms_elem = metabolite.find("hmdb:synonyms", self.namespace)
            logp = self.get_experimental_properties(metabolite, "logp")
            state = self.get_text(metabolite, "state")

            subject_node = {
                "id": primary_id,
                "name": name.lower() if name else None,
                "synonym": self.get_list(synonyms_elem, "synonym")
                if synonyms_elem is not None
                else [],
                "description": description,
                "chemical_formula": self.get_text(metabolite, "chemical_formula"),
                "molecular_weight": self.get_molecular_weights(metabolite),
                "state": state.lower() if state else None,
                "water_solubility": self.get_experimental_properties(
                    metabolite, "water_solubility"
                ),
                "logp": logp,
                "melting_point": self.get_experimental_properties(metabolite, "melting_point"),
                "type": "biolink:SmallMolecule",
                "xrefs": xrefs,
            }
            subject_node = self.remove_empty_none_values(subject_node)

            diseases_elem = metabolite.find("hmdb:diseases", self.namespace)
            if diseases_elem is not None:
                for disease_elem in diseases_elem.findall("hmdb:disease", self.namespace):
                    references = self.get_medi_references(disease_elem)
                    if references is None:
                        continue

                    association_node = {
                        "predicate": "biolink:ChemicalEntityToDiseaseAssociation",
                        "infores": "hmdb_v5",
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


if __name__ == "__main__":
    zip_path = os.path.join("downloads", "hmdb_metabolites.zip")
    hmdb_xml = extract_xml_from_zip(zip_path)
    parser = HMDBParse(hmdb_xml)
    
    mime_records = [record for record in parser.parse_mime()]
    save_pickle(mime_records, "hmdb_v5_microbe_metabolite.pkl")
    for record in mime_records:
        print(record[0])

    medi_records = [record for record in parser.parse_medi()]
    save_pickle(medi_records, "hmdb_v5_metabolite_disease.pkl")
    for record in medi_records:
        print(record[0])
