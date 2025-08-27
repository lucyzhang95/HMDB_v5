"""Record management for HMDB data parsing and output generation."""

import json
import os
import sys
import time
from itertools import chain
from pathlib import Path
from typing import Dict, Iterator, List

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'parsers'))

from tqdm.auto import tqdm

from .cache_helper import load_pickle, save_json, save_pickle
from .cache_manager import CacheManager
from metabolite_parser import HMDBMetaboliteParser
from protein_parser import HMDBProteinParser
from .reader import extract_file_from_zip


class RecordManager:
    """Manages the generation and storage of HMDB association records."""

    def __init__(self, data_dir: str = "downloads"):
        self.data_dir = Path(data_dir).resolve()

    def _prepare_xml_files(self) -> tuple[str, str]:
        """Prepare XML files by extracting from ZIP if necessary."""
        # handle HMDB metabolites XML
        metabolite_xml = self.data_dir / "hmdb_metabolites.xml"
        if not metabolite_xml.is_file():
            metabolite_zip = self.data_dir / "hmdb_metabolites.zip"
            if not metabolite_zip.exists():
                raise FileNotFoundError(f"Neither {metabolite_xml} nor {metabolite_zip} found")
            metabolite_xml = extract_file_from_zip(str(metabolite_zip), "hmdb_metabolites.xml")

        # handle HMDB proteins XML
        protein_xml = self.data_dir / "hmdb_proteins.xml"
        if not protein_xml.is_file():
            protein_zip = self.data_dir / "hmdb_proteins.zip"
            if not protein_zip.exists():
                raise FileNotFoundError(f"‼️Neither {protein_xml} nor {protein_zip} found")
            protein_xml = extract_file_from_zip(str(protein_zip), "hmdb_proteins.xml")

        return str(metabolite_xml), str(protein_xml)

    def generate_all_records(self) -> Dict[str, List[Dict]]:
        """Generate all HMDB association records."""
        print("⚙️ Generating HMDB association records...")

        # step 1: prepare XML files
        metabolite_xml, protein_xml = self._prepare_xml_files()

        # step 2: initialize parsers
        metabolite_parser = HMDBMetaboliteParser(metabolite_xml)
        protein_parser = HMDBProteinParser(protein_xml)

        # step 3: generate records with progress tracking
        records = {}

        print("▶️ Parsing microbe-metabolite associations...")
        records["microbe-metabolite"] = list(
            tqdm(metabolite_parser.parse_microbe_metabolite(), desc="Microbe-Metabolite")
        )

        print("▶️ Parsing metabolite-disease associations...")
        records["metabolite-disease"] = list(
            tqdm(metabolite_parser.parse_metabolite_disease(), desc="Metabolite-Disease")
        )

        print("▶️ Parsing metabolite-protein associations...")
        records["metabolite-protein"] = list(
            tqdm(metabolite_parser.parse_metabolite_protein(), desc="Metabolite-Protein")
        )

        print("▶️ Parsing metabolite-pathway associations...")
        records["metabolite-pathway"] = list(
            tqdm(metabolite_parser.parse_metabolite_pathway(), desc="Metabolite-Pathway")
        )

        print("▶️ Parsing protein-pathway associations...")
        records["protein-pathway"] = list(
            tqdm(protein_parser.parse_protein_pathway(), desc="Protein-Pathway")
        )

        print("▶️ Parsing protein-biological process associations...")
        records["protein-biological_process"] = list(
            tqdm(
                protein_parser.parse_protein_biological_process(), desc="Protein-BiologicalProcess"
            )
        )

        return records

    def save_records(self, records: Dict[str, List[Dict]], output_format: str = "both") -> Dict:
        """Save records to cache files."""
        print("Saving parsed records...")

        if output_format in ("pickle", "both"):
            save_pickle(records, "hmdb_v5_parsed_records.pkl")

        if output_format in ("json", "both"):
            save_json(records, "hmdb_v5_parsed_records.json")

        stats = self._calculate_record_stats(records)
        save_json(stats, "hmdb_v5_record_stats.json")

        return stats

    def _calculate_record_stats(self, records: Dict[str, List[Dict]]) -> Dict:
        """Calculate statistics for generated records."""
        stats = {
            "total_records": sum(len(record_list) for record_list in records.values()),
            "association_types": {},
            "generation_timestamp": time.time(),
        }

        for assoc_type, record_list in records.items():
            stats["association_types"][assoc_type] = {
                "count": len(record_list),
                "percentage": len(record_list) / stats["total_records"] * 100
                if stats["total_records"] > 0
                else 0,
            }

        return stats

    def load_cached_records(self) -> Dict[str, List[Dict]]:
        """Load cached records if they exist."""
        cached = load_pickle("hmdb_v5_parsed_records.pkl")
        if cached:
            return cached
        else:
            raise FileNotFoundError("‼️No cached records found. Run --generate-records first.")

    def get_record_iterator(self) -> Iterator[Dict]:
        """Get an iterator over all records."""
        try:
            records = self.load_cached_records()
            return chain(*records.values())
        except FileNotFoundError:
            print("❌ No cached records found. Generating records...")
            records = self.generate_all_records()
            self.save_records(records)
            return chain(*records.values())

    def export_records(
        self, output_file: str, format: str = "json", association_types: List[str] = None
    ) -> None:
        """Export records to a specified format."""
        records = self.load_cached_records()

        if association_types:
            filtered_records = {k: v for k, v in records.items() if k in association_types}
        else:
            filtered_records = records

        output_path = Path("cache") / Path(output_file).name
        output_path.parent.mkdir(parents=True, exist_ok=True)

        if format.lower() == "jsonl":
            self._export_jsonl(filtered_records, output_path)
        elif format.lower() == "json":
            self._export_json(filtered_records, output_path)
        elif format.lower() == "tsv":
            self._export_tsv(filtered_records, output_path)
        else:
            raise ValueError(f"Unsupported format: {format}")

        print(f"Records exported to {output_path}")

    def _export_jsonl(self, records: Dict, output_path: Path) -> None:
        """Export records as JSON Lines format."""
        with open(output_path, "w") as f:
            for record_list in records.values():
                for record in record_list:
                    f.write(json.dumps(record) + "\n")

    def _export_json(self, records: Dict, output_path: Path) -> None:
        """Export records as JSON format."""
        with open(output_path, "w") as f:
            json.dump(records, f, indent=2)

    def _export_tsv(self, records: Dict, output_path: Path) -> None:
        """Export records as TSV format (flattened)."""
        import csv

        all_fields = set()
        all_records = list(chain(*records.values()))

        for record in all_records:
            all_fields.update(self._flatten_dict(record).keys())

        fieldnames = sorted(all_fields)

        with open(output_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            for record in all_records:
                flattened = self._flatten_dict(record)
                writer.writerow(flattened)

    def _flatten_dict(self, d: Dict, parent_key: str = "", sep: str = ".") -> Dict:
        """Flatten nested dictionary for TSV export."""
        items = []
        for k, v in d.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            if isinstance(v, dict):
                items.extend(self._flatten_dict(v, new_key, sep=sep).items())
            elif isinstance(v, list):
                items.append((new_key, json.dumps(v)))
            else:
                items.append((new_key, v))
        return dict(items)

    def get_record_summary(self) -> Dict:
        """Get summary statistics of cached records."""
        try:
            records = self.load_cached_records()
            stats_file = "hmdb_v5_record_stats.json"

            if os.path.exists(os.path.join("cache", stats_file)):
                with open(os.path.join("cache", stats_file)) as f:
                    stats = json.load(f)
            else:
                stats = self._calculate_record_stats(records)
                save_json(stats, stats_file)

            return stats

        except FileNotFoundError:
            return {"error": "No cached records found", "total_records": 0, "association_types": {}}


def cache_hmdb_database(
    email: str, umls_api_key: str, data_dir: str = "downloads", force_refresh: bool = False
) -> Dict:
    """
    Complete HMDB database caching and record generation pipeline.

    :param email: Email for NCBI Entrez API
    :param umls_api_key: API key for UMLS
    :param data_dir: Directory containing HMDB XML files
    :param force_refresh: Whether to force refresh of all cached data
    
    :return: Dictionary with pipeline results and statistics
    """
    start_time = time.time()

    print("HMDB Database Caching Pipeline")
    print("=" * 50)

    cache_manager = CacheManager(email, umls_api_key)

    if force_refresh or not cache_manager.is_cache_complete():
        print("Caching reference data (taxonomies, diseases, proteins)...")

        data_path = Path(data_dir).resolve()

        metabolite_xml = data_path / "hmdb_metabolites.xml"
        if not metabolite_xml.is_file():
            metabolite_zip = data_path / "hmdb_metabolites.zip"
            metabolite_xml = extract_file_from_zip(str(metabolite_zip), "hmdb_metabolites.xml")

        protein_xml = data_path / "hmdb_proteins.xml"
        if not protein_xml.is_file():
            protein_zip = data_path / "hmdb_proteins.zip"
            protein_xml = extract_file_from_zip(str(protein_zip), "hmdb_proteins.xml")

        cache_manager.cache_metabolite_data(str(metabolite_xml))
        cache_manager.cache_protein_data(str(protein_xml))
    else:
        print("Reference data cache is complete, skipping...")

    print("\nGenerating association records...")
    record_manager = RecordManager(data_dir)

    if force_refresh or not os.path.exists(os.path.join("cache", "hmdb_v5_parsed_records.pkl")):
        records = record_manager.generate_all_records()
        stats = record_manager.save_records(records)
    else:
        print("Parsed records already exist, loading from cache...")
        stats = record_manager.get_record_summary()

    duration = time.time() - start_time

    print("\n" + "=" * 50)
    print("HMDB Database Caching Complete!")
    print(f"Total time: {duration / 60:.2f} minutes")
    print(f"Total association records: {stats.get('total_records', 0):,}")

    if "association_types" in stats:
        print("\nAssociation breakdown:")
        for assoc_type, type_stats in stats["association_types"].items():
            print(f"  {assoc_type}: {type_stats['count']:,} ({type_stats['percentage']:.1f}%)")

    return {
        "success": True,
        "duration": duration,
        "total_records": stats.get("total_records", 0),
        "association_types": stats.get("association_types", {}),
        "cache_complete": cache_manager.is_cache_complete(),
    }
