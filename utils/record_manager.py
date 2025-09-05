"""Record management for HMDB data parsing and output generation."""

import json
import logging
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any, Dict, List, Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "parsers"))

from metabolite_parser import HMDBMetaboliteParser
from protein_parser import HMDBProteinParser
from tqdm.auto import tqdm

from .cache_manager import CacheManager
from .reader import extract_file_from_zip
from .record_deduplicator import _create_fingerprint, deduplicate_and_merge

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class RecordManager:
    """Manages the generation and storage of HMDB association records."""

    def __init__(self, data_dir: str = "downloads"):
        """
        Initialize RecordManager.

        :param data_dir: Directory containing HMDB data raw files
        """
        self.data_dir = Path(data_dir).resolve()
        if not self.data_dir.exists():
            raise FileNotFoundError(f"Data directory does not exist: {self.data_dir}")

    def _prepare_xml_files(self) -> Tuple[str, str]:
        """
        Prepare XML files by extracting from ZIP if necessary.

        :return: Tuple of (metabolite_xml_path, protein_xml_path)
        :raises FileNotFoundError: If no XML files are found in the data directory
        """
        metabolite_xml = self._get_xml_file("hmdb_metabolites")
        protein_xml = self._get_xml_file("hmdb_proteins")

        return str(metabolite_xml), str(protein_xml)

    def _get_xml_file(self, base_name: str) -> Path:
        """
        Get XML file path, extracting from ZIP if necessary.

        :param base_name: Base name without extension (e.g., "hmdb_metabolites")
        :return: Path to the XML file
        :raises FileNotFoundError: If no XML files are found in the data directory
        """
        xml_path = self.data_dir / f"{base_name}.xml"

        if xml_path.is_file():
            return xml_path

        zip_path = self.data_dir / f"{base_name}.zip"
        if not zip_path.exists():
            raise FileNotFoundError(
                f"Neither {xml_path} nor {zip_path} found. "
                f"Please ensure HMDB data files are downloaded to {self.data_dir}"
            )

        logger.info(f"Extracting {base_name}.xml from {zip_path}")
        extracted_path = extract_file_from_zip(str(zip_path), f"{base_name}.xml")
        return Path(extracted_path)

    def _flatten_dict(
            self, d: Dict[str, Any], parent_key: str = "", sep: str = "."
    ) -> Dict[str, Any]:
        """
        Flatten nested dictionary for JSONL export.

        :param d: Dictionary to flatten
        :param parent_key: the Parent key for nested structure
        :param sep: Separator for nested keys
        :return: Flattened dictionary
        """
        items = []
        for k, v in d.items():
            new_key = f"{parent_key}{sep}{k}" if parent_key else k
            if isinstance(v, dict):
                items.extend(self._flatten_dict(v, new_key, sep=sep).items())
            elif isinstance(v, list):
                items.append((new_key, json.dumps(v) if v else ""))
            else:
                items.append((new_key, str(v) if v is not None else ""))
        return dict(items)

    def _get_parser_tasks(
            self, metabolite_parser: HMDBMetaboliteParser, protein_parser: HMDBProteinParser
    ) -> List[Tuple[str, Any, str]]:
        """
        Get a list of parser tasks with their configurations.

        :param metabolite_parser: Metabolite parser instance
        :param protein_parser: Protein parser instance
        
        :return: List of (key, parser_function, description) tuples
        """
        return [
            (
                "microbe-metabolite",
                metabolite_parser.parse_microbe_metabolite,
                "Microbe-Metabolite",
            ),
            (
                "metabolite-disease",
                metabolite_parser.parse_metabolite_disease,
                "Metabolite-Disease",
            ),
            (
                "metabolite-protein",
                metabolite_parser.parse_metabolite_protein,
                "Metabolite-Protein",
            ),
            (
                "metabolite-pathway",
                metabolite_parser.parse_metabolite_pathway,
                "Metabolite-Pathway",
            ),
            ("protein-pathway", protein_parser.parse_protein_pathway, "Protein-Pathway"),
            (
                "protein-biological_process",
                protein_parser.parse_protein_biological_process,
                "Protein-BiologicalProcess",
            ),
        ]

    def _get_task_configs(self, tasks: List[Tuple[str, Any, str]]) -> List[Tuple[str, Any, str]]:
        """Return task configurations for lazy iterator creation."""
        return tasks

    def generate_and_export_streamed(self, output_file: str, output_format: str = "jsonl") -> None:
        """
        Generate, deduplicate, and export records in a memory-efficient way.

        :param output_file: Path to the output file
        :param output_format: Output format (json, jsonl, pkl)
        :raises ValueError: If output_format is not supported
        :raises IOError: If file operations fail
        """
        SUPPORTED_FORMATS = {"json", "jsonl", "pkl"}
        BATCH_SIZE = 1000

        if output_format not in SUPPORTED_FORMATS:
            raise ValueError(
                f"!! Unsupported format '{output_format}'. Must be one of: {SUPPORTED_FORMATS}"
            )

        overall_start_time = time.time()
        logger.info(
            f"\n>>> Starting generation and export to '{output_file}' (format: {output_format})"
        )

        output_path = Path("records") / Path(output_file).name
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # initialize parsers
        try:
            setup_start = time.time()
            metabolite_xml, protein_xml = self._prepare_xml_files()
            metabolite_parser = HMDBMetaboliteParser(metabolite_xml)
            protein_parser = HMDBProteinParser(protein_xml)
            setup_time = time.time() - setup_start
            logger.info(f"-> Parser setup completed in {setup_time:.2f} seconds")
        except Exception as e:
            logger.error(f"!!! Failed to prepare parsers: {e}")
            raise

        tasks = self._get_parser_tasks(metabolite_parser, protein_parser)

        # initialize statistics
        stats = {
            "total_records": 0,
            "association_types": {},
            "generation_timestamp": time.time(),
            "processing_times": {"setup": setup_time},
        }
        raw_counts = {}

        # process records
        try:
            processed_association_types = self._process_and_export_records(
                output_path, output_format, tasks, stats, raw_counts, BATCH_SIZE
            )

            if output_format == "json":
                self._combine_json_files(output_path, processed_association_types)

        except Exception as e:
            logger.error(f"!!! Failed to process records: {e}")
            raise

        # finalize statistics
        self._finalize_statistics(stats, raw_counts)

        # print summary
        overall_time = time.time() - overall_start_time
        stats["processing_times"]["total_pipeline"] = overall_time
        self._print_export_summary(raw_counts, stats, output_path, overall_time)

    def _process_and_export_records(
            self,
            output_path: Path,
            output_format: str,
            tasks: list,
            stats: dict,
            raw_counts: dict,
            batch_size: int,
    ) -> List[str]:
        """Process and export records with proper resource management."""

        task_configs = self._get_task_configs(tasks)
        processed_association_types = []

        with tqdm(
                task_configs, desc="Processing association types", unit="type", position=0, leave=True
        ) as overall_pbar:

            for task_index, (key, parser_func, desc) in enumerate(overall_pbar):
                try:
                    logger.info(f">>> Creating iterator for {desc}...")
                    record_iterator = parser_func()

                    if output_format == "json":
                        temp_output_path = output_path.parent / f"{key}.json"
                    elif output_format == "pkl":
                        temp_output_path = output_path.parent / f"{key}.pkl"
                    else:
                        temp_output_path = output_path

                    total_raw_count, total_deduped_count = self._process_single_association_type(
                        temp_output_path,
                        key,
                        record_iterator,
                        desc,
                        output_format,
                        batch_size,
                        overall_pbar,
                        is_first_association_type=(task_index == 0),
                    )

                    if total_deduped_count > 0:
                        processed_association_types.append(key)
                        raw_counts[key] = total_raw_count
                        stats["total_records"] += total_deduped_count
                        stats["association_types"][key] = {"count": total_deduped_count}

                except Exception as e:
                    logger.error(f"!!! Error processing {desc}: {e}")
                    continue

        return processed_association_types

    def _process_single_association_type(
            self,
            output_path: Path,
            key: str,
            record_iterator,
            desc: str,
            output_format: str,
            batch_size: int,
            overall_pbar,
            is_first_association_type: bool = True,
    ) -> Tuple[int, int]:
        """Process a single association type with streaming and batched processing."""
        association_start = time.time()
        overall_pbar.set_description(f"Processing {desc}")
        logger.info(f">>> Processing: {desc}")

        # process records with batched deduplication
        total_raw_count, total_deduped_count = self._process_records_streamed(
            output_path, record_iterator, key, output_format, batch_size, is_first_association_type
        )

        processing_time = time.time() - association_start

        if total_deduped_count == 0:
            logger.info(f"!! Skipping empty association type: {key}")
            return 0, 0

        # log completion
        duplicates_removed = total_raw_count - total_deduped_count
        logger.info(
            f"[DONE] {desc} completed: "
            f"{total_deduped_count:,} records ({duplicates_removed:,} duplicates removed) "
            f"in {processing_time:.2f}s (streamed processing)"
        )

        return total_raw_count, total_deduped_count

    def _process_records_streamed(
            self,
            output_path: Path,
            record_iterator,
            key: str,
            output_format: str,
            batch_size: int,
            is_first_association_type: bool = True,
    ) -> Tuple[int, int]:
        """
        Process records in true streaming fashion with batched deduplication.

        This processes records in small batches to minimize memory usage while
        still allowing for deduplication within reasonable memory constraints.

        Returns:
            Tuple of (total_raw_count, total_deduped_count)
        """
        total_raw_count = 0
        total_deduped_count = 0
        batch = []
        global_seen_keys = set()  # track duplicates across batches

        # initialize progress tracking
        record_pbar = tqdm(desc=f"Processing {key}", unit="records", position=1, leave=False)

        if output_format == "pkl":
            file_context = open(output_path, "wb")
        else:
            if output_format == "jsonl" and not is_first_association_type:
                file_mode = "a"
            else:
                file_mode = "w"
            file_context = open(output_path, file_mode, newline="", encoding="utf-8")

        try:
            with file_context as f:
                if output_format == "json" and is_first_association_type:
                    f.write("[\n")

                first_record_written = False

                for record in record_iterator:
                    batch.append(record)
                    total_raw_count += 1
                    record_pbar.update(1)
                    record_pbar.set_description(f"Processing {key} ({total_raw_count:,} records)")

                    if len(batch) >= batch_size:
                        deduped_count = self._process_batch(
                            f, batch, key, output_format, global_seen_keys, first_record_written
                        )
                        total_deduped_count += deduped_count
                        if deduped_count > 0:
                            first_record_written = True
                        batch.clear()

                if batch:
                    deduped_count = self._process_batch(
                        f, batch, key, output_format, global_seen_keys, first_record_written
                    )
                    total_deduped_count += deduped_count
                    batch.clear()

                if output_format == "json":
                    f.write("\n]")

        finally:
            record_pbar.close()

        return total_raw_count, total_deduped_count

    def _process_batch(
            self,
            f,
            batch: list,
            key: str,
            output_format: str,
            global_seen_keys: set,
            first_record_written: bool,
    ) -> int:
        """
        Process a single batch: deduplicate and write records.

        Returns:
            Number of records written from this batch
        """
        if not batch:
            return 0

        # deduplicate within a batch
        deduped_batch = deduplicate_and_merge(batch)

        new_records_no_dup = []
        for record in deduped_batch:
            record_fingerprint = _create_fingerprint(record)
            if record_fingerprint not in global_seen_keys:
                global_seen_keys.add(record_fingerprint)
                new_records_no_dup.append(record)

        if new_records_no_dup:
            if output_format == "pkl":
                for record in new_records_no_dup:
                    pickle.dump(record, f)
            elif output_format == "jsonl":
                for record in new_records_no_dup:
                    record["association_type"] = key
                    json.dump(record, f, ensure_ascii=False)
                    f.write("\n")
            elif output_format == "json":
                for i, record in enumerate(new_records_no_dup):
                    record["association_type"] = key

                    if first_record_written or i > 0:
                        f.write(",\n")

                    json_str = json.dumps(record, ensure_ascii=False, indent=2)
                    indented_lines = []
                    for line in json_str.split("\n"):
                        indented_lines.append("  " + line)
                    f.write("\n".join(indented_lines))

                    first_record_written = True

        return len(new_records_no_dup)

    def _combine_json_files(self, output_path: Path, association_types: List[str]) -> None:
        """
        Combine individual JSON files into the final grouped JSON structure.

        :param output_path: Path to the output JSON file.
        :param association_types: List of association type keys that were processed
        """
        logger.info("\n>>> Combining JSON files into final structure...")

        temp_final_path = output_path.with_suffix(".tmp")

        with open(temp_final_path, "w", encoding="utf-8") as final_file:
            final_file.write("{\n")

            written_types = 0
            for assoc_type in association_types:
                temp_file = output_path.parent / f"{assoc_type}.json"

                if temp_file.exists() and temp_file.stat().st_size > 2:
                    if written_types > 0:
                        final_file.write(",\n")

                    final_file.write(f'  "{assoc_type}": ')

                    with open(temp_file, "r", encoding="utf-8") as temp:
                        final_file.write(temp.read())

                    written_types += 1

                # free temp file
                if temp_file.exists():
                    temp_file.unlink()

            final_file.write("\n}")

        temp_final_path.replace(output_path)
        logger.info(f"[DONE] Combined {written_types} association types into {output_path}")

    def _finalize_statistics(self, stats: dict, raw_counts: dict) -> None:
        """Calculate final statistics including percentages."""
        total = stats["total_records"]
        for key in stats["association_types"]:
            count = stats["association_types"][key]["count"]
            stats["association_types"][key]["percentage"] = (
                (count / total * 100) if total > 0 else 0
            )

    def _print_export_summary(
            self,
            raw_counts: Dict[str, int],
            stats: Dict[str, Any],
            output_path: Path,
            overall_time: float,
    ) -> None:
        """Print export summary statistics with detailed timing information."""
        logger.info("[DONE] Deduplication and export complete!")

        logger.info("\nProcessing Summary:")
        logger.info("-" * 80)
        for key in raw_counts:
            dedup_count = stats.get("association_types", {}).get(key, {}).get("count", 0)
            removed = raw_counts[key] - dedup_count

            logger.info(f"  {key:25} | " f"Records: {dedup_count:>8,} | " f"Removed: {removed:>8,}")

        # print timing breakdown
        logger.info("-" * 80)
        logger.info(f"-> Total pipeline:  {overall_time:>8.2f}s")

        total_records = stats.get("total_records", 0)
        if overall_time > 0:
            records_per_sec = total_records / overall_time
            logger.info(f"-> Processing rate: {records_per_sec:>8,.1f} records/second")

        logger.info("-" * 80)
        logger.info(f"\nSuccess! {total_records:,} records exported to {output_path}")
        logger.info(f"-> File size: {self._get_file_size(output_path)}")

    def _get_file_size(self, file_path: Path) -> str:
        """Get human-readable file size."""
        try:
            size_bytes = file_path.stat().st_size
            if size_bytes == 0:
                return "0 B"

            size_names = ["B", "KB", "MB", "GB", "TB"]
            import math

            i = int(math.floor(math.log(size_bytes, 1024)))
            p = math.pow(1024, i)
            s = round(size_bytes / p, 2)
            return f"{s} {size_names[i]}"
        except OSError:
            return "Unknown size"


def cache_hmdb_database(
        email: str, umls_api_key: str, data_dir: str = "downloads", force_refresh: bool = False
) -> Dict[str, Any]:
    """
    HMDB reference data caching pipeline (does not generate association records).

    This function only caches reference data (taxonomies, diseases, proteins, etc.)
    which is required before generating association records. Use the separate
    generate_and_export_streamed() method to create association records.

    Args:
        email: Email for NCBI Entrez API
        umls_api_key: API key for UMLS
        data_dir: Directory containing HMDB XML files
        force_refresh: Whether to force refresh of all cached data

    Returns:
        Dictionary with pipeline results (no record statistics)
    """
    start_time = time.time()

    logger.info("\nHMDB Reference Data Caching Pipeline")
    logger.info("=" * 50)

    try:
        cache_manager = CacheManager(email, umls_api_key)

        if force_refresh or not cache_manager.is_cache_complete():
            logger.info(">>> Caching reference data (taxonomies, diseases, proteins)...")

            record_manager = RecordManager(data_dir)
            metabolite_xml, protein_xml = record_manager._prepare_xml_files()

            # cache reference data with progress tracking
            metabolite_start = time.time()
            logger.info("️>>> Caching metabolite reference data...")
            cache_manager.cache_metabolite_data(metabolite_xml)
            metabolite_time = time.time() - metabolite_start
            logger.info(f"[DONE] Metabolite data cached in {metabolite_time:.2f} seconds")

            protein_start = time.time()
            logger.info(">>> Caching protein reference data...")
            cache_manager.cache_protein_data(protein_xml)
            protein_time = time.time() - protein_start
            logger.info(f"[DONE] Protein data cached in {protein_time:.2f} seconds")

        else:
            logger.info("[DONE] Reference data cache is complete, skipping...")
            metabolite_time = 0
            protein_time = 0

        duration = time.time() - start_time

        logger.info("=" * 50)
        logger.info("[DONE] HMDB Reference Data Caching Complete!")
        logger.info(f"-> Total time: {duration:.2f} seconds ({duration / 60:.2f} minutes)")

        if force_refresh or not cache_manager.is_cache_complete():
            logger.info(
                f"-> Breakdown: Metabolite: {metabolite_time:.2f}s, Protein: {protein_time:.2f}s"
            )

        logger.info("\n Next steps:")
        logger.info("   Use generate_and_export_streamed() to create association records")
        logger.info(
            "   Example: record_manager.generate_and_export_streamed('output.jsonl', 'jsonl')"
        )

        return {
            "success": True,
            "duration": duration,
            "cache_complete": cache_manager.is_cache_complete(),
            "operations_performed": {
                "metabolite_caching": force_refresh or not cache_manager.is_cache_complete(),
                "protein_caching": force_refresh or not cache_manager.is_cache_complete(),
            },
            "timing": {
                "metabolite_caching": metabolite_time,
                "protein_caching": protein_time,
                "total": duration,
            },
        }

    except Exception as e:
        logger.error(f"!!! Reference data caching pipeline failed: {e}")
        return {
            "success": False,
            "error": str(e),
            "duration": time.time() - start_time,
            "cache_complete": False,
        }
