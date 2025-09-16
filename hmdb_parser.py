"""
HMDB Parser for v5 - Modular HMDB data extraction and processing pipeline.

This is the main entry point for the restructured HMDB parser that extracts
biomedical associations from HMDB metabolite and protein XML files.
"""

import argparse
import logging
import os
import sys
from pathlib import Path
from typing import Dict, Iterator, Optional, Tuple

sys.path.append(os.path.join(os.path.dirname(__file__), "utils"))

from utils.cache_pipeline import show_cache_status
from utils.record_manager import RecordManager, cache_hmdb_database

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def load_hmdb_data(data_dir: str = "downloads") -> Iterator[Dict]:
    """
    Load HMDB association data with automatic caching.

    Args:
        data_dir: Directory containing HMDB data files

    Yields:
        Dictionary records of HMDB associations

    Raises:
        FileNotFoundError: If data directory or required files don't exist
    """
    try:
        record_manager = RecordManager(data_dir)
        yield from record_manager.get_record_iterator()
    except Exception as e:
        logger.error(f"Failed to load HMDB data: {e}")
        raise


def validate_environment_variables(
        email: Optional[str], umls_api_key: Optional[str]
) -> Tuple[str, str]:
    """Validate required environment variables."""
    if not email:
        raise ValueError(
            "Email is required. Set EMAIL_ADDRESS environment variable or use --email argument."
        )
    if not umls_api_key:
        raise ValueError(
            "UMLS API key is required. Set UMLS_API_KEY environment variable or use --umls-api-key argument."
        )
    return email, umls_api_key


def setup_argument_parser():
    """Set up and return the argument parser."""

    parser = argparse.ArgumentParser(
        description="HMDB Parser - Extract biomedical associations from HMDB data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Usage Examples:
    # Show current cache status
    python hmdb_parser.py --status

    # Cache reference data, then generate and export records to JSONL (recommended workflow)
    python hmdb_parser.py --cache-reference-data
    python hmdb_parser.py --export hmdb_v5_parsed_records.jsonl --format jsonl

    # Do everything in one go (cache reference data and then export)
    python hmdb_parser.py --cache-reference-data --export hmdb_v5_parsed_records.jsonl --format jsonl

    # Export only specific formats
    python hmdb_parser.py --export hmdb_v5_parsed_records.json --format json

    # Force refresh all cached data
    python hmdb_parser.py --cache-all --force-refresh
        """,
    )

    parser.add_argument(
        "--email",
        default=os.getenv("EMAIL_ADDRESS"),
        help="Email for NCBI Entrez API (can also use EMAIL_ADDRESS env var)",
    )
    parser.add_argument(
        "--umls-api-key",
        default=os.getenv("UMLS_API_KEY"),
        help="UMLS API key (can also use UMLS_API_KEY env var)",
    )

    parser.add_argument(
        "--data-dir",
        default="downloads",
        help="Directory containing HMDB XML files (default: downloads)",
    )

    cache_group = parser.add_argument_group("caching operations")
    cache_group.add_argument(
        "--cache-all", action="store_true", help="Run complete caching pipeline for reference data"
    )
    cache_group.add_argument(
        "--cache-reference-data",
        action="store_true",
        help="Cache only reference data (taxonomies, diseases, etc.)",
    )
    cache_group.add_argument(
        "--status", action="store_true", help="Show cache status and statistics"
    )

    export_group = parser.add_argument_group("export operations")
    export_group.add_argument("--export", help="Generate and export records to the specified file")
    export_group.add_argument(
        "--format",
        choices=["json", "jsonl", "pkl"],
        default="json",
        help="Export format (default: jsonl, recommended for large datasets)",
    )
    export_group.add_argument(
        "--association-types",
        nargs="*",
        help="Filter export to specific association types (not supported in streaming mode)",
    )

    parser.add_argument(
        "--force-refresh", action="store_true", help="Force refresh of all cached data"
    )
    parser.add_argument("--quiet", action="store_true", help="Reduce output verbosity")
    parser.add_argument(
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Set logging level (default: INFO)",
    )

    return parser


def configure_logging(log_level: str, quiet: bool) -> None:
    """Configure logging based on arguments."""
    if quiet:
        level = logging.WARNING
    else:
        level = getattr(logging, log_level.upper())

    logging.getLogger().setLevel(level)
    logging.getLogger(__name__).setLevel(level)


def main():
    """Main command-line interface."""
    from dotenv import load_dotenv

    load_dotenv()

    parser = setup_argument_parser()
    args = parser.parse_args()

    configure_logging(args.log_level, args.quiet)

    try:
        email, umls_api_key = validate_environment_variables(args.email, args.umls_api_key)
    except ValueError as e:
        logger.error(str(e))
        parser.error(str(e))

    data_path = Path(args.data_dir)
    if not data_path.exists():
        logger.error(f"!!! Data directory does not exist: {data_path}")
        parser.error(f"!!! Data directory does not exist: {data_path}")

    try:
        if args.status:
            show_cache_status(email, umls_api_key, args.data_dir)

        if args.cache_all or args.cache_reference_data:
            logger.info("Running HMDB reference data caching pipeline...")
            result = cache_hmdb_database(
                email=email,
                umls_api_key=umls_api_key,
                data_dir=args.data_dir,
                force_refresh=args.force_refresh,
            )

            if result.get("success"):
                logger.info("Reference data caching completed successfully.")
                if not args.quiet:
                    logger.info(f"Duration: {result.get('duration', 0) / 60:.2f} minutes")
            else:
                logger.error(f"Caching failed: {result.get('error', 'Unknown error')}")
                sys.exit(1)

        if args.export:
            logger.info(f"Exporting records to {args.export} in {args.format} format...")
            record_manager = RecordManager(args.data_dir)

            try:
                record_manager.generate_and_export_streamed(
                    output_file=args.export, output_format=args.format
                )
                logger.info("[DONE] Export completed successfully.")
            except Exception as e:
                logger.error(f"!!! Export failed: {e}")
                if not args.quiet:
                    import traceback

                    traceback.print_exc()
                sys.exit(1)

        elif not any([args.status, args.cache_all, args.cache_reference_data, args.export]):
            logger.info("No action specified. Use --help for usage information.")
            parser.print_help()

    except KeyboardInterrupt:
        logger.info("Operation cancelled by user.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"An unexpected error occurred: {e}")
        if not args.quiet:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
    # use "python hmdb_parser.py --export hmdb_v5_parsed_records.jsonl --format jsonl"
