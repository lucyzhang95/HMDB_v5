"""Manual taxon name to taxonomy ID mappings for HMDB data."""

# manual corrections for text2term mappings
TEXT2TERM_CORRECTIONS = {
    "geobacillus thermoglucosidasius": 1426,
    "bacteroides spp.": 816,
    "clostridium difficile": 1496,
    "pseudomonas pseudomaleii": 28450,
    "lactobacillus plantarum": 1590,
    "clostridium butylicum": 1492,
    "corynebacterium jekeium": 38289,
}

# exclusions from text2term mappings
TEXT2TERM_EXCLUSIONS = {
    "gram-negative bacteria",
}

# manual taxon name to taxonomy ID mappings
MANUAL_TAXON_MAPPINGS = {
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
    "biï¿dobacterium": 1678,
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


def get_text2term_corrections() -> dict:
    """Get text2term mapping corrections."""
    return TEXT2TERM_CORRECTIONS.copy()


def get_text2term_exclusions() -> set:
    """Get items to exclude from text2term mappings."""
    return TEXT2TERM_EXCLUSIONS.copy()


def get_manual_taxon_mappings() -> dict:
    """Get manual taxon name to taxonomy ID mappings."""
    return MANUAL_TAXON_MAPPINGS.copy()


def apply_text2term_corrections(text2term_mapped: dict) -> dict:
    """Apply manual corrections to text2term mappings."""
    corrections = get_text2term_corrections()
    exclusions = get_text2term_exclusions()

    # apply corrections
    for name, corrected_taxid in corrections.items():
        if name in text2term_mapped:
            text2term_mapped[name]["taxid"] = corrected_taxid
            text2term_mapped[name]["mapping_tool"] = "manual"

    # remove exclusions
    for exclusion in exclusions:
        text2term_mapped.pop(exclusion, None)

    return text2term_mapped


def format_manual_mappings(taxon_names: list[str]) -> dict:
    """Format manual mappings for given taxon names."""
    manual_mappings = get_manual_taxon_mappings()
    return {
        name: {"taxid": manual_mappings[name], "mapping_tool": "manual"}
        for name in set(taxon_names)
        if name in manual_mappings
    }
