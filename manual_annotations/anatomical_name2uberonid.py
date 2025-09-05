"""Module for correcting specific anatomical names to their corresponding UBERON IDs."""

UBERON_CORRECTIONS = {
    "all tissues": {
        "id": "UBERON:0000479",
        "name": "tissue",
        "original_name": "all tissues",
        "category": "biolink:AnatomicalEntity",
    },
    "smooth muscle": {
        "id": "UBERON:0001135",
        "name": "smooth muscle tissue",
        "original_name": "smooth muscle",
        "category": "biolink:AnatomicalEntity",
    },
}

MANUAL_ANATOMY_MAPPINGS = {
    "neuron": {
        "id": "CL:0000540",
        "name": "neuron",
        "original_name": "neuron",
        "category": "biolink:Cell",
    },
    "erythrocyte": {
        "id": "CL:0000232",
        "name": "erythrocyte",
        "original_name": "erythrocyte",
        "category": "biolink:Cell",
    },
    "leukocyte": {
        "id": "CL:0000738",
        "name": "leukocyte",
        "original_name": "leukocyte",
        "category": "biolink:Cell",
    },
    "platelet": {
        "id": "CL:0000233",
        "name": "platelet",
        "original_name": "platelet",
        "category": "biolink:Cell",
    },
    "cellular cytoplasm": {
        "id": "GO:0005737",
        "name": "cytoplasm",
        "original_name": "cellular cytoplasm",
        "category": "biolink:CellularComponent",
    },
    "pericardial effusion": {
        "id": "UBERON:0002406",
        "name": "pericardial sac",
        "original_name": "pericardial effusion",
        "category": "biolink:AnatomicalEntity",
    },
}


def _get_uberon_correction():
    return UBERON_CORRECTIONS.copy()


def _get_uberon_manual_mappings():
    return MANUAL_ANATOMY_MAPPINGS.copy()


def _apply_uberon_corrections(uberon_mapped: dict) -> dict:
    correction = _get_uberon_correction()
    for name, corrected_data in correction.items():
        if name in uberon_mapped:
            uberon_mapped[name] = corrected_data
    return uberon_mapped


def apply_manual_anatomy_mappings(uberon_mapped: dict) -> dict:
    uberon_corrected = _apply_uberon_corrections(uberon_mapped)
    uberon_corrected.update(MANUAL_ANATOMY_MAPPINGS)
    return uberon_corrected
