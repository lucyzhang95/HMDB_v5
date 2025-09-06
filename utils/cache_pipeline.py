"""High-level cache pipeline orchestration for HMDB data processing."""

import time
from pathlib import Path

from .cache_helper import clear_cache, get_cache_size
from .cache_manager import CacheManager
from .reader import extract_file_from_zip


class CachePipeline:
    """High-level orchestration of HMDB data caching pipeline."""

    def __init__(self, email: str, umls_api_key: str, data_dir: str = "downloads"):
        self.email = email
        self.umls_api_key = umls_api_key
        self.data_dir = Path(data_dir).resolve()
        self.cache_manager = CacheManager(email, umls_api_key)

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
                raise FileNotFoundError(f"Neither {protein_xml} nor {protein_zip} found")
            protein_xml = extract_file_from_zip(str(protein_zip), "hmdb_proteins.xml")

        return str(metabolite_xml), str(protein_xml)

    def run_full_cache_pipeline(
            self, force_refresh: bool = False, skip_existing: bool = True
    ) -> dict:
        """Run the complete caching pipeline."""
        start_time = time.time()

        print("HMDB Data Caching Pipeline")
        print("=" * 50)

        if force_refresh:
            print("Force refresh requested - clearing existing cache...")
            clear_cache()

        try:
            metabolite_xml, protein_xml = self._prepare_xml_files()
        except FileNotFoundError as e:
            print(f"Error: {e}")
            return {"success": False, "error": str(e)}

        if skip_existing and self.cache_manager.is_cache_complete():
            print("[DONE] All cache files already exist. Skipping caching process.")
            print("Use force_refresh=True to rebuild cache.")
            return {
                "success": True,
                "message": "Cache already complete",
                "duration": time.time() - start_time,
            }

        missing_files = self.cache_manager.get_missing_cache_files()
        if missing_files:
            print(f"!!! Missing cache files: {len(missing_files)}")
            for file in missing_files[:5]:
                print(f"  - {file}")
            if len(missing_files) > 5:
                print(f"  ... and {len(missing_files) - 5} more")

        try:
            # cache metabolite data
            if not self.cache_manager.cache_exists_for_file(metabolite_xml) or force_refresh:
                self.cache_manager.cache_metabolite_data(metabolite_xml)
            else:
                print("[DONE] Metabolite cache exists, skipping...")

            # cache protein data
            if not self.cache_manager.cache_exists_for_file(protein_xml) or force_refresh:
                self.cache_manager.cache_protein_data(protein_xml)
            else:
                print("[DONE] Protein cache exists, skipping...")

            duration = time.time() - start_time

            print("\n" + "=" * 50)
            print("[DONE] Cache Pipeline Complete!")
            print(f"-> Total time: {duration / 60:.2f} minutes")

            cache_info = get_cache_size()
            total_size_mb = sum(info["size_mb"] for info in cache_info.values())
            print(f"-> Total cache size: {total_size_mb:.1f} MB")
            print(f"-> Number of cache files: {len(cache_info)}")

            return {
                "success": True,
                "duration": duration,
                "cache_size_mb": total_size_mb,
                "num_files": len(cache_info),
            }

        except Exception as e:
            print(f"\n!!! Error during caching: {str(e)}")
            return {"success": False, "error": str(e)}

    def validate_cache_integrity(self) -> dict:
        """Validate the integrity of cached data."""
        print("\n>>> Validating cache integrity...")

        status = self.cache_manager.get_cache_status()
        missing = [f for f, exists in status.items() if not exists]

        validation_results = {
            "files_missing": missing,
            "total_files": len(status),
            "files_present": len(status) - len(missing),
        }

        if not missing:
            try:
                cached_data = self.cache_manager.load_cached_data()
                validation_results["data_loadable"] = True
                validation_results["taxon_count"] = len(cached_data.get("taxon_info", {}))
                validation_results["disease_count"] = len(cached_data.get("disease_info", {}))
                validation_results["pathway_count"] = len(
                    cached_data.get("pathway_descriptions", {})
                )
            except Exception as e:
                validation_results["data_loadable"] = False
                validation_results["load_error"] = str(e)

        return validation_results

    def show_cache_summary(self) -> None:
        """Display a summary of the current cache state."""
        print("\nHMDB Cache Summary")
        print("=" * 30)

        status = self.cache_manager.get_cache_status()
        cache_info = get_cache_size()

        print(f"-> Cache files: {sum(status.values())}/{len(status)} present")

        if cache_info:
            total_size = sum(info["size_mb"] for info in cache_info.values())
            print(f"-> Total size: {total_size:.1f} MB")

            largest_files = sorted(cache_info.items(), key=lambda x: x[1]["size_mb"], reverse=True)[:5]

            print("\nLargest cache files:")
            for filename, info in largest_files:
                print(f"  {filename}: {info['size_mb']:.1f} MB")

        missing = [f for f, exists in status.items() if not exists]
        if missing:
            print(f"\n!!! Missing files ({len(missing)}):")
            for filename in missing[:5]:
                print(f"  - {filename}")
            if len(missing) > 5:
                print(f"  ... and {len(missing) - 5} more")

        # validate data if complete
        if not missing:
            print("\n>>> Validating cache integrity...")
            try:
                cached_data = self.cache_manager.load_cached_data()
                print("[DONE] All cache files load successfully")
                print(f"  - Taxon mappings: {len(cached_data.get('taxon_info', {}))}")
                print(f"  - Disease mappings: {len(cached_data.get('disease_info', {}))}")
                print(
                    f"  - Pathway descriptions: {len(cached_data.get('pathway_descriptions', {}))}"
                )
            except Exception as e:
                print(f"!!! Cache validation failed: {e}")


def run_cache_pipeline(
        email: str, umls_api_key: str, data_dir: str = "downloads", force_refresh: bool = False
) -> dict:
    """Convenience function to run the complete cache pipeline."""
    pipeline = CachePipeline(email, umls_api_key, data_dir)
    return pipeline.run_full_cache_pipeline(force_refresh=force_refresh)


def show_cache_status(email: str, umls_api_key: str, data_dir: str = "downloads") -> None:
    """Convenience function to show cache status."""
    pipeline = CachePipeline(email, umls_api_key, data_dir)
    pipeline.show_cache_summary()
