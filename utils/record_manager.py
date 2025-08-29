"""Record management for HMDB data parsing and output generation."""

import csv
import json
import logging
import os
import pickle
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), ".."))
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "parsers"))

import warnings

from metabolite_parser import HMDBMetaboliteParser
from protein_parser import HMDBProteinParser
from tqdm.auto import tqdm

from .cache_manager import CacheManager
from .reader import extract_file_from_zip
from .record_deduplicator import deduplicate_and_merge

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
        Flatten nested dictionary for TSV export.

        :param d: Dictionary to flatten
        :param parent_key: Parent key for nested structure
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
        Get list of parser tasks with their configurations.

        Args:
            metabolite_parser: Metabolite parser instance
            protein_parser: Protein parser instance

        Returns:
            List of (key, parser_function, description) tuples
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

    def _scan_headers_for_tsv(self, tasks: List[Tuple[str, Any, str]]) -> List[str]:
        """
        Scan all records to determine TSV headers.

        Args:
            tasks: List of parser tasks

        Returns:
            Sorted list of unique field names
        """
        logger.info("Scanning all records to determine TSV headers (this may take a while)...")
        all_headers = set()

        for _, parser_func, desc in tqdm(tasks, desc="Header Scan"):
            try:
                for record in parser_func():
                    record["association_type"] = desc
                    all_headers.update(self._flatten_dict(record).keys())
            except Exception as e:
                logger.warning(f"Error scanning headers for {desc}: {e}")
                continue

        fieldnames = sorted(list(all_headers))
        logger.info(f"Found {len(fieldnames)} unique headers")
        return fieldnames

    def _get_task_iterators(
        self, tasks: List[Tuple[str, Any, str]]
    ) -> List[Tuple[str, Iterator, str]]:
        """
        Get task iterators without pre-scanning (truly streaming approach).

        Args:
            tasks: List of parser tasks

        Returns:
            List of tasks with iterators (empty tasks will be filtered during processing)
        """
        task_iterators = []

        for key, parser_func, desc in tasks:
            try:
                iterator = parser_func()
                task_iterators.append((key, iterator, desc))
            except Exception as e:
                logger.error(f"Error creating iterator for {key}: {e}")
                continue

        return task_iterators

    def _write_json_record(
        self, file_handle, record: Dict[str, Any], is_first: bool, is_last: bool
    ) -> None:
        """Write a single record in JSON format."""
        json.dump(record, file_handle, indent=2)
        if not is_last:
            file_handle.write(",\n")

    def _write_records_by_format(
        self,
        file_handle,
        records: List[Dict[str, Any]],
        key: str,
        output_format: str,
        writer: Optional[csv.DictWriter] = None,
        is_first_type: bool = True,
        processed_count: int = 0,
    ) -> None:
        """
        Write records in the specified format with progress tracking.

        Args:
            file_handle: Output file handle
            records: List of records to write
            key: Association type key
            output_format: Output format (json, jsonl, tsv)
            writer: CSV writer for TSV format
            is_first_type: Whether this is the first association type
            processed_count: Number of association types already processed
        """
        record_iterator = tqdm(
            enumerate(records),
            total=len(records),
            desc=f"Exporting {key}",
            unit="records",
            leave=False,
        )

        for i, record in record_iterator:
            record["association_type"] = key

            if output_format == "jsonl":
                file_handle.write(json.dumps(record) + "\n")
            elif output_format == "tsv" and writer:
                writer.writerow(self._flatten_dict(record))
            elif output_format == "json":
                if processed_count == 0 and i == 0:
                    file_handle.write(f'  "{key}": [\n')
                elif i == 0:
                    file_handle.write(f',\n  "{key}": [\n')

                self._write_json_record(file_handle, record, i == 0, i == len(records) - 1)

                if i == len(records) - 1:
                    file_handle.write("\n  ]")

    def generate_and_export_streamed(self, output_file: str, output_format: str = "jsonl") -> None:
        """
        Generate, deduplicate, and export records in a memory-efficient way.

        Args:
            output_file: Output file path
            output_format: Export format (json, jsonl, pkl, tsv)
        """
        overall_start_time = time.time()
        logger.info(
            f"🚀 Starting generation and export to '{output_file}' (format: {output_format})"
        )

        if output_format in ("pkl", "tsv"):
            warnings.warn(
                f"The '{output_format}' format can be slow or memory-intensive. "
                "'jsonl' is recommended for large datasets.",
                UserWarning,
            )

        # Ensure cache directory exists
        output_path = Path("cache") / Path(output_file).name
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            setup_start = time.time()
            metabolite_xml, protein_xml = self._prepare_xml_files()
            metabolite_parser = HMDBMetaboliteParser(metabolite_xml)
            protein_parser = HMDBProteinParser(protein_xml)
            setup_time = time.time() - setup_start
            logger.info(f"⚙️ Parser setup completed in {setup_time:.2f} seconds")
        except Exception as e:
            logger.error(f"Failed to prepare parsers: {e}")
            raise

        tasks = self._get_parser_tasks(metabolite_parser, protein_parser)

        # Initialize statistics
        stats = {
            "total_records": 0,
            "association_types": {},
            "generation_timestamp": time.time(),
            "processing_times": {},
        }
        raw_counts = {}
        all_records_for_pickle = {}

        # Handle TSV header scanning if needed
        fieldnames = []
        if output_format == "tsv":
            header_start = time.time()
            fieldnames = self._scan_headers_for_tsv(tasks)
            header_time = time.time() - header_start
            logger.info(f"📋 Header scanning completed in {header_time:.2f} seconds")
            stats["processing_times"]["header_scan"] = header_time

        try:
            with open(output_path, "w", newline="", encoding="utf-8") as f:
                writer = None
                if output_format == "tsv":
                    writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
                    writer.writeheader()
                elif output_format == "json":
                    f.write("{\n")

                task_iterators = self._get_task_iterators(tasks)

                # Create overall progress bar for all association types
                overall_pbar = tqdm(
                    task_iterators,
                    desc="Processing association types",
                    unit="type",
                    position=0,
                    leave=True,
                )

                processed_tasks = 0  # Counter for valid tasks (for JSON formatting)

                for _, (key, record_iterator, desc) in enumerate(overall_pbar):
                    association_start = time.time()
                    overall_pbar.set_description(f"Processing {desc}")

                    try:
                        # Parse records with progress bar - collect to list to check if empty
                        logger.info(f"📊 Processing: {desc}")
                        parsing_start = time.time()

                        # Process records in batches to check if empty without loading everything
                        raw_records = []
                        batch_size = 1000
                        record_count = 0

                        # Create a tqdm iterator that we can break early if empty
                        record_pbar = tqdm(
                            record_iterator,
                            desc=f"Parsing {key}",
                            unit="records",
                            position=1,
                            leave=False,
                        )

                        for record in record_pbar:
                            raw_records.append(record)
                            record_count += 1

                            # Update progress bar with current count
                            record_pbar.set_description(f"Parsing {key} ({record_count:,} records)")

                            # Process in batches to manage memory
                            if len(raw_records) >= batch_size:
                                # Continue processing - we have data
                                pass

                        record_pbar.close()

                        # Skip if no records found
                        if not raw_records:
                            logger.info(f"⚠️  Skipping empty association type: {key}")
                            continue

                        parsing_time = time.time() - parsing_start
                        raw_counts[key] = len(raw_records)

                        # Deduplicate with progress tracking
                        dedup_start = time.time()
                        logger.info(f"🔄 Deduplicating {len(raw_records):,} records for {key}...")
                        deduped_records = deduplicate_and_merge(raw_records)
                        dedup_time = time.time() - dedup_start
                        count = len(deduped_records)

                        # Export records
                        export_start = time.time()
                        if output_format == "pkl":
                            all_records_for_pickle[key] = deduped_records
                        else:
                            self._write_records_by_format(
                                f,
                                deduped_records,
                                key,
                                output_format,
                                writer,
                                is_first_type=(processed_tasks == 0),
                                processed_count=processed_tasks,
                            )
                        export_time = time.time() - export_start

                        # Update statistics
                        association_total_time = time.time() - association_start
                        stats["total_records"] += count
                        stats["association_types"][key] = {"count": count}
                        stats["processing_times"][key] = {
                            "parsing": parsing_time,
                            "deduplication": dedup_time,
                            "export": export_time,
                            "total": association_total_time,
                        }

                        processed_tasks += 1

                        # Log timing info
                        duplicates_removed = raw_counts[key] - count
                        logger.info(
                            f"✅ {desc} completed: "
                            f"{count:,} records ({duplicates_removed:,} duplicates removed) "
                            f"in {association_total_time:.2f}s "
                            f"(parse: {parsing_time:.2f}s, dedup: {dedup_time:.2f}s, export: {export_time:.2f}s)"
                        )

                        # Clean up memory
                        del raw_records
                        if output_format != "pkl":
                            del deduped_records

                    except Exception as e:
                        logger.error(f"❌ Error processing {desc}: {e}")
                        continue

                overall_pbar.close()

                if output_format == "json":
                    f.write("\n}\n")

        except IOError as e:
            logger.error(f"Failed to write output file: {e}")
            raise

        # save pickle format separately
        if output_format == "pkl":
            pickle_start = time.time()
            logger.info("💾 Saving all collected records to pickle file...")
            with open(output_path, "wb") as pkl_f:
                pickle.dump(all_records_for_pickle, pkl_f)
            pickle_time = time.time() - pickle_start
            stats["processing_times"]["pickle_save"] = pickle_time
            logger.info(f"💾 Pickle save completed in {pickle_time:.2f} seconds")

        # update statistics
        total = stats["total_records"]
        for key in stats["association_types"]:
            count = stats["association_types"][key]["count"]
            stats["association_types"][key]["percentage"] = (
                (count / total * 100) if total > 0 else 0
            )

        # Print summary
        overall_time = time.time() - overall_start_time
        stats["processing_times"]["total_pipeline"] = overall_time
        self._print_export_summary(raw_counts, stats, output_path, overall_time)

    def _print_export_summary(
        self,
        raw_counts: Dict[str, int],
        stats: Dict[str, Any],
        output_path: Path,
        overall_time: float,
    ) -> None:
        """Print export summary statistics with detailed timing information."""
        logger.info("🎉 Deduplication and export complete!")

        # Print per-association type summary
        logger.info("\n📊 Processing Summary:")
        logger.info("-" * 80)
        for key in raw_counts:
            dedup_count = stats.get("association_types", {}).get(key, {}).get("count", 0)
            removed = raw_counts[key] - dedup_count

            timing = stats.get("processing_times", {}).get(key, {})
            parse_time = timing.get("parsing", 0)
            dedup_time = timing.get("deduplication", 0)
            export_time = timing.get("export", 0)
            total_time = timing.get("total", 0)

            logger.info(
                f"  {key:25} | "
                f"Records: {dedup_count:>8,} | "
                f"Removed: {removed:>8,} | "
                f"Time: {total_time:>6.2f}s "
                f"(P:{parse_time:>5.2f}s D:{dedup_time:>5.2f}s E:{export_time:>5.2f}s)"
            )

        # Print timing breakdown
        logger.info("-" * 80)
        processing_times = stats.get("processing_times", {})

        if "header_scan" in processing_times:
            logger.info(f"⏱️  Header scanning: {processing_times['header_scan']:>8.2f}s")
        if "pickle_save" in processing_times:
            logger.info(f"⏱️  Pickle saving:   {processing_times['pickle_save']:>8.2f}s")

        total_processing = sum(
            timing.get("total", 0)
            for timing in processing_times.values()
            if isinstance(timing, dict)
        )
        logger.info(f"⏱️  Data processing: {total_processing:>8.2f}s")
        logger.info(f"⏱️  Total pipeline:  {overall_time:>8.2f}s")

        # Performance metrics
        total_records = stats.get("total_records", 0)
        if overall_time > 0:
            records_per_sec = total_records / overall_time
            logger.info(f"⚡ Processing rate: {records_per_sec:>8,.1f} records/second")

        logger.info("-" * 80)
        logger.info(f"✅ Success! {total_records:,} records exported to {output_path}")
        logger.info(f"📁 File size: {self._get_file_size(output_path)}")

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

    logger.info("🗃️  HMDB Reference Data Caching Pipeline")
    logger.info("=" * 50)

    try:
        cache_manager = CacheManager(email, umls_api_key)

        if force_refresh or not cache_manager.is_cache_complete():
            logger.info("📚 Caching reference data (taxonomies, diseases, proteins)...")

            # Prepare XML files
            record_manager = RecordManager(data_dir)
            metabolite_xml, protein_xml = record_manager._prepare_xml_files()

            # Cache reference data with progress tracking
            metabolite_start = time.time()
            logger.info("🧬 Caching metabolite reference data...")
            cache_manager.cache_metabolite_data(metabolite_xml)
            metabolite_time = time.time() - metabolite_start
            logger.info(f"✅ Metabolite data cached in {metabolite_time:.2f} seconds")

            protein_start = time.time()
            logger.info("🧪 Caching protein reference data...")
            cache_manager.cache_protein_data(protein_xml)
            protein_time = time.time() - protein_start
            logger.info(f"✅ Protein data cached in {protein_time:.2f} seconds")

        else:
            logger.info("✅ Reference data cache is complete, skipping...")
            metabolite_time = 0
            protein_time = 0

        duration = time.time() - start_time

        logger.info("=" * 50)
        logger.info("🎉 HMDB Reference Data Caching Complete!")
        logger.info(f"⏱️  Total time: {duration:.2f} seconds ({duration / 60:.2f} minutes)")

        if force_refresh or not cache_manager.is_cache_complete():
            logger.info(
                f"📊 Breakdown: Metabolite: {metabolite_time:.2f}s, Protein: {protein_time:.2f}s"
            )

        logger.info("\n💡 Next steps:")
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
        logger.error(f"❌ Reference data caching pipeline failed: {e}")
        return {
            "success": False,
            "error": str(e),
            "duration": time.time() - start_time,
            "cache_complete": False,
        }
