"""XML reading utilities for HMDB data extraction."""

import os
import pathlib
import zipfile
from typing import Iterator, Union

from lxml import etree as ET


def extract_file_from_zip(zip_path: str, expected_filename: str) -> str:
    """Extract a file from a ZIP archive if it doesn't already exist."""
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
    """Strip namespace from XML element tag."""
    return tag.split("}", 1)[-1] if "}" in tag else tag


def get_all_microbe_names(input_xml: Union[str, pathlib.Path]) -> Iterator[str]:
    """Extract microbe taxon names from HMDB metabolites XML."""
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
    """Extract all diseases and their OMIM IDs from HMDB metabolites XML."""
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


def get_all_uniprot_ids_from_hmdb(input_xml: Union[str, pathlib.Path]) -> dict[str, str]:
    """Extract protein name to UniProt ID mapping from HMDB metabolites XML."""
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    protein2uniprot = {}
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
                protein2uniprot[protein_name] = uniprot_id if uniprot_id else None
    return protein2uniprot


def get_all_uniprot_ids_from_hmdbp(input_xml: Union[str, pathlib.Path]) -> list[str]:
    """Extract all UniProt IDs from HMDB proteins XML."""
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


def get_all_go_terms_from_hmdbp(input_xml: Union[str, pathlib.Path]) -> list[str]:
    """Extract all GO terms from HMDB proteins XML."""
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


def get_all_anatomical_terms_from_hmdb(input_xml: Union[str, pathlib.Path]) -> list[str]:
    """Extract all tissue terms from HMDB metabolites XML."""
    namespace = {"hmdb": "http://www.hmdb.ca"}
    tree = ET.parse(input_xml)
    root = tree.getroot()

    anatomical_terms = set()

    for metabolite in root.findall("hmdb:metabolite", namespace):
        bio_prop = metabolite.find("hmdb:biological_properties", namespace)
        if bio_prop is not None:
            tissue_elem = bio_prop.find("hmdb:tissue_locations", namespace)
            if tissue_elem is not None:
                for tissue in tissue_elem.findall("hmdb:tissue", namespace):
                    if tissue.text and tissue.text.strip():
                        anatomical_terms.add(tissue.text.strip().lower())

            specimen_elem = bio_prop.find("hmdb:biospecimen_locations", namespace)
            if specimen_elem is not None:
                for specimen in specimen_elem.findall("hmdb:biospecimen", namespace):
                    if specimen.text and specimen.text.strip():
                        anatomical_terms.add(specimen.text.strip().lower())

    return sorted(list(anatomical_terms))
