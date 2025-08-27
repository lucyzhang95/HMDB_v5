"""
HMDB Parser for v5 - Modular HMDB data extraction and processing pipeline.

This is the main entry point for the restructured HMDB parser that extracts
biomedical associations from HMDB metabolite and protein XML files.
"""

import os
import sys
from typing import Dict, Iterator

sys.path.append(os.path.join(os.path.dirname(__file__), "utils"))


from utils.cache_pipeline import CachePipeline, show_cache_status
from utils.record_manager import RecordManager, cache_hmdb_database


def load_hmdb_data(data_dir: str = "downloads") -> Iterator[Dict]:
    """
    Load HMDB association data with automatic caching.

    This is the main function equivalent to the original load_hmdb_data(),
    but now uses the modular architecture with proper caching.

    :param data_dir: Directory containing HMDB XML files

    :yield: HMDB association records in standardized format
    """
    record_manager = RecordManager(data_dir)
    yield from record_manager.get_record_iterator()


def main():
    """Main command-line interface."""
    import argparse

    from dotenv import load_dotenv

    load_dotenv()

    parser = argparse.ArgumentParser(
        description="HMDB v5 Parser - Extract biomedical associations from HMDB data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Usage Examples:
    # Cache all reference data and generate association records
    python hmdb_parser.py --cache-all

    # Show current cache status
    python hmdb_parser.py --status

    # Force refresh of all cached data
    python hmdb_parser.py --cache-all --force-refresh

    # Generate records only (requires cached reference data)
    python hmdb_parser.py --generate-records

    # Export records to different formats
    python hmdb_parser.py --export records.jsonl --format jsonl
    python hmdb_parser.py --export records.json --format json
    python hmdb_parser.py --export records.tsv --format tsv
    
    # Example Workflow: 1. cache reference data, 2. generate records, and 3. export
    python hmdb_parser.py --cache-reference-data --generate-records --export records.json --format json
    
    or
    
    # python hmdb_parser.py --cache-all --export records.json --format json
        """,
    )

    parser.add_argument(
        "--email",
        default=os.getenv("EMAIL_ADDRESS"),
        help="Email for NCBI Entrez API (or set EMAIL_ADDRESS env var)",
    )
    parser.add_argument(
        "--umls-api-key",
        default=os.getenv("UMLS_API_KEY"),
        help="UMLS API key (or set UMLS_API_KEY env var)",
    )

    parser.add_argument(
        "--data-dir",
        default="downloads",
        help="Directory containing HMDB XML files (default: downloads)",
    )

    parser.add_argument(
        "--cache-all",
        action="store_true",
        help="Run complete caching pipeline (reference data + records)",
    )
    parser.add_argument(
        "--cache-reference-data",
        action="store_true",
        help="Cache only reference data (taxonomies, diseases, proteins)",
    )
    parser.add_argument(
        "--generate-records",
        action="store_true",
        help="Generate association records (requires cached reference data)",
    )
    parser.add_argument("--status", action="store_true", help="Show cache status and statistics")

    parser.add_argument("--export", help="Export records to specified file")
    parser.add_argument(
        "--format",
        choices=["json", "jsonl", "tsv"],
        default="jsonl",
        help="Export format (default: jsonl)",
    )
    parser.add_argument(
        "--association-types", nargs="*", help="Filter export to specific association types"
    )

    parser.add_argument(
        "--force-refresh", action="store_true", help="Force refresh of all cached data"
    )
    parser.add_argument("--quiet", action="store_true", help="Reduce output verbosity")

    args = parser.parse_args()

    if not args.email:
        print("Error: Email is required. Use --email or set EMAIL_ADDRESS environment variable.")
        sys.exit(1)

    if not args.umls_api_key:
        print(
            "Error: UMLS API key is required. Use --umls-api-key or set UMLS_API_KEY environment variable."
        )
        sys.exit(1)

    try:
        if args.status:
            show_cache_status(args.email, args.umls_api_key, args.data_dir)

        elif args.cache_all:
            print("Running complete HMDB caching pipeline...")
            result = cache_hmdb_database(
                email=args.email,
                umls_api_key=args.umls_api_key,
                data_dir=args.data_dir,
                force_refresh=args.force_refresh,
            )

            if result["success"]:
                print(f"\nPipeline completed successfully in {result['duration'] / 60:.1f} minutes")
                print(f"Generated {result['total_records']:,} association records")
            else:
                print(f"Pipeline failed: {result.get('error', 'Unknown error')}")
                sys.exit(1)

        elif args.cache_reference_data:
            print("Caching reference data only...")
            pipeline = CachePipeline(args.email, args.umls_api_key, args.data_dir)

            # prepare XML files
            metabolite_xml, protein_xml = pipeline._prepare_xml_files()

            # cache reference data
            pipeline.cache_manager.cache_metabolite_data(metabolite_xml)
            pipeline.cache_manager.cache_protein_data(protein_xml)

            print("Reference data caching complete!")

        elif args.generate_records:
            print("Generating association records...")
            record_manager = RecordManager(args.data_dir)
            records = record_manager.generate_all_records()
            stats = record_manager.save_records(records)

            print(f"Generated {stats['total_records']:,} records")
            for assoc_type, type_stats in stats["association_types"].items():
                print(f"  {assoc_type}: {type_stats['count']:,}")

        elif args.export:
            print(f"Exporting records to {args.export}...")
            record_manager = RecordManager(args.data_dir)
            record_manager.export_records(
                args.export, format=args.format, association_types=args.association_types
            )

        else:
            print("No action specified. Use --help for usage information.")
            parser.print_help()

    except Exception as e:
        print(f"Error: {e}")
        if not args.quiet:
            import traceback

            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
