"""Module for correcting specific anatomical names to their corresponding UBERON IDs."""

OXO_CORRECTIONS = {
    "all tissues": {
        "id": "UBERON:0000479",
        "name": "tissue",
        "original_name": "all tissues",
        "type": "biolink:AnatomicalEntity",
    },
    "smooth muscle": {
        "id": "UBERON:0001135",
        "name": "smooth muscle tissue",
        "original_name": "smooth muscle",
        "type": "biolink:AnatomicalEntity",
    },
}


def get_oxo_correction():
    return OXO_CORRECTIONS.copy()


def apply_oxo_correction(uberon_mapped: dict) -> dict:
    correction = get_oxo_correction()
    for name, corrected_data in correction.items():
        if name in uberon_mapped:
            uberon_mapped[name] = corrected_data
    return uberon_mapped
