"""
Manual annotations and corrections for HMDB data mapping.

This module contains manually curated mappings and corrections for
taxonomy and disease name resolution that couldn't be automatically
resolved through standard ontology services.
"""

from .disease_name2id import format_manual_disease_mappings, get_manual_disease_mappings
from .taxon_name2taxid import (
    apply_text2term_corrections,
    format_manual_mappings,
    get_manual_taxon_mappings,
    get_text2term_corrections,
    get_text2term_exclusions,
)

__all__ = [
    "get_text2term_corrections",
    "get_text2term_exclusions",
    "get_manual_taxon_mappings",
    "apply_text2term_corrections",
    "format_manual_mappings",
    "get_manual_disease_mappings",
    "format_manual_disease_mappings",
]
