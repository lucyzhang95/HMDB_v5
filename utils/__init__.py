"""
HMDB Parser Utils - Modular components for HMDB data processing.
"""

import os
import sys

from .cache_helper import (
    cache_exists,
    clear_cache,
    get_cache_path,
    get_cache_size,
    list_cache_files,
    load_pickle,
    save_json,
    save_pickle,
)
from .cache_manager import CacheManager
from .cache_pipeline import CachePipeline, run_cache_pipeline, show_cache_status
from .ontology_mapper import DiseaseMapper, ProteinMapper, TaxonMapper
from .ontology_services import (
    BiothingsServices,
    DiseaseServices,
    GeneServices,
    GOServices,
    NCITServices,
    PathwayServices,
    ProteinServices,
    TaxonServices,
    UMLSClient,
)
from .parser_helper import (
    IDHierarchy,
    MolecularWeightExtractor,
    OrganismClassifier,
    PropertyExtractor,
    ReferenceExtractor,
    XMLParseHelper,
)

# core modules
from .reader import (
    extract_file_from_zip,
    get_all_diseases,
    get_all_go_terms_from_hmdbp,
    get_all_microbe_names,
    get_all_uniprot_ids_from_hmdb,
    get_all_uniprot_ids_from_hmdbp,
)
from .record_manager import RecordManager, cache_hmdb_database

sys.path.append(os.path.join(os.path.dirname(__file__), "..", "parsers"))

from metabolite_parser import HMDBMetaboliteParser
from protein_parser import HMDBProteinParser

__version__ = "1.0"
__author__ = "Lucy Zhang (BioThings Team)"
__description__ = "Modular HMDB data extraction and processing pipeline"


# convenience imports for backward compatibility
def load_hmdb_data(data_dir="downloads"):
    """Load HMDB data with automatic caching (backward compatibility)."""
    from .record_manager import RecordManager

    record_manager = RecordManager(data_dir)
    yield from record_manager.get_record_iterator()


__all__ = [
    # reader functions
    "extract_file_from_zip",
    "get_all_microbe_names",
    "get_all_diseases",
    "get_all_uniprot_ids_from_hmdb",
    "get_all_uniprot_ids_from_hmdbp",
    "get_all_go_terms_from_hmdbp",
    # cache functions
    "save_pickle",
    "load_pickle",
    "save_json",
    "cache_exists",
    "get_cache_path",
    "clear_cache",
    "get_cache_size",
    "list_cache_files",
    # cache management
    "CacheManager",
    "CachePipeline",
    "run_cache_pipeline",
    "show_cache_status",
    # ontology mapping
    "TaxonMapper",
    "DiseaseMapper",
    "ProteinMapper",
    # services
    "TaxonServices",
    "DiseaseServices",
    "ProteinServices",
    "GeneServices",
    "GOServices",
    "UMLSClient",
    "NCITServices",
    "BiothingsServices",
    "PathwayServices",
    # parser helpers
    "XMLParseHelper",
    "OrganismClassifier",
    "ReferenceExtractor",
    "IDHierarchy",
    "MolecularWeightExtractor",
    "PropertyExtractor",
    # parsers
    "HMDBMetaboliteParser",
    "HMDBProteinParser",
    # record management
    "RecordManager",
    "cache_hmdb_database",
    "load_hmdb_data",
]
