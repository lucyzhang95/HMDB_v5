import glob
import json
import os
import random
from collections import Counter, defaultdict
from datetime import datetime
from typing import Any, Dict, List, Optional, Set


class HMDBRecordStatsReporter:
    """Generates comprehensive statistics for HMDB parsed records."""

    def __init__(self):
        self.report_dir = os.path.join("..", "reports")
        os.makedirs(self.report_dir, exist_ok=True)
        print("HMDBRecordStatsReporter initialized.")

        self.all_subject_properties = set()
        self.all_object_properties = set()
        self.all_association_properties = set()

    def _find_latest_jsonl_file(self, base_path: str = "hmdb_v5_parsed_records") -> Optional[str]:
        """Find the most recent timestamped JSONL file."""
        pattern = os.path.join("..", "records", f"{base_path}_*.jsonl")
        files = glob.glob(pattern)

        if not files:
            print(f"No JSONL files found matching pattern: {pattern}")
            return None

        latest_file = max(files, key=os.path.getmtime)
        print(f"Found latest JSONL file: {latest_file}")
        return latest_file

    def _load_jsonl_records(self, file_path: str):
        """Generator to load records one by one from JSONL file."""
        print(f"Loading records from: {file_path}")

        record_count = 0
        with open(file_path, "r", encoding="utf-8") as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    try:
                        record = json.loads(line)
                        record_count += 1
                        yield record
                    except json.JSONDecodeError as e:
                        print(f"Error parsing line {line_num}: {e}")
                        continue

        print(f"Total records processed: {record_count}")

    def _extract_all_properties(self, data: Dict[str, Any], prefix: str = "") -> Set[str]:
        """Recursively extract all property paths from a dictionary."""
        properties = set()

        if not isinstance(data, dict):
            return properties

        for key, value in data.items():
            current_path = f"{prefix}.{key}" if prefix else key
            properties.add(current_path)

            if isinstance(value, dict):
                nested_props = self._extract_all_properties(value, current_path)
                properties.update(nested_props)
            elif isinstance(value, list) and value and isinstance(value[0], dict):
                nested_props = self._extract_all_properties(value[0], current_path)
                properties.update(nested_props)

        return properties

    def _count_property_stats(
            self, records: List[Dict[str, Any]], property_path: str
    ) -> Dict[str, Any]:
        """Count statistics for a specific property path."""
        values = []
        non_null_count = 0

        for record in records:
            value = self._get_nested_value(record, property_path)
            if value is not None:
                non_null_count += 1
                if isinstance(value, list):
                    values.extend([str(v) for v in value if v is not None])
                else:
                    values.append(str(value))

        unique_values = set(values)
        total_records = len(records)

        return {
            "raw_count": len(values),
            "record_count": non_null_count,
            "unique_count": len(unique_values),
            "percentage": round((non_null_count / total_records) * 100, 2)
            if total_records > 0
            else 0,
            "sample_values": list(unique_values)[:10],  # First 10 unique values as sample
        }

    def _get_nested_value(self, data: Dict[str, Any], path: str) -> Any:
        """Get value from nested dictionary using dot notation path."""
        keys = path.split(".")
        current = data

        for key in keys:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return None

        return current

    def _analyze_duplicates(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze duplicate records by _id and export samples."""
        print("Analyzing duplicate records...")

        id_to_records = defaultdict(list)
        for record in records:
            record_id = record.get("_id")
            if record_id:
                id_to_records[record_id].append(record)

        duplicate_groups = {}
        duplicate_count_distribution = Counter()

        for record_id, record_list in id_to_records.items():
            if len(record_list) > 1:
                duplicate_groups[record_id] = record_list
                duplicate_count_distribution[len(record_list)] += 1

        sampled_duplicates = {}
        for count, _ in duplicate_count_distribution.items():
            ids_with_count = [rid for rid, rlist in duplicate_groups.items() if len(rlist) == count]

            sample_size = min(3, len(ids_with_count))
            sampled_ids = random.sample(ids_with_count, sample_size)

            sampled_duplicates[count] = {}
            for record_id in sampled_ids:
                sampled_duplicates[count][record_id] = duplicate_groups[record_id]

        self._export_duplicate_records(duplicate_groups, sampled_duplicates)

        duplicate_stats = {
            "total_duplicate_groups": len(duplicate_groups),
            "total_duplicate_records": sum(len(rlist) for rlist in duplicate_groups.values()),
            "duplication_count_distribution": dict(duplicate_count_distribution),
            "max_duplicates_for_single_id": max(duplicate_count_distribution.keys())
            if duplicate_count_distribution
            else 0,
            "sampled_duplicates_by_count": {
                count: list(group_data.keys()) for count, group_data in sampled_duplicates.items()
            },
        }

        return duplicate_stats

    def _export_duplicate_records(
            self,
            all_duplicates: Dict[str, List[Dict]],
            sampled_duplicates: Dict[int, Dict[str, List[Dict]]],
    ) -> str:
        """Export duplicate records to JSON files."""
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        sampled_export_path = os.path.join(
            self.report_dir, f"sampled_duplicate_records_{timestamp}.json"
        )
        sampled_export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "description": "Randomly sampled duplicate records (max 3 per duplication count)",
                "total_sampled_groups": sum(len(groups) for groups in sampled_duplicates.values()),
            },
            "duplicates_by_count": sampled_duplicates,
        }

        with open(sampled_export_path, "w", encoding="utf-8") as f:
            json.dump(sampled_export_data, f, indent=2, sort_keys=True, default=str)

        all_export_path = os.path.join(self.report_dir, f"{timestamp}_all_duplicate_records.json")
        all_export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "description": "All duplicate records found",
                "total_duplicate_groups": len(all_duplicates),
                "total_duplicate_records": sum(len(rlist) for rlist in all_duplicates.values()),
            },
            "duplicate_records": all_duplicates,
        }

        with open(all_export_path, "w", encoding="utf-8") as f:
            json.dump(all_export_data, f, indent=2, sort_keys=True, default=str)

        print(f"Sampled duplicate records exported to: {sampled_export_path}")
        print(f"All duplicate records exported to: {all_export_path}")

        return sampled_export_path

    def _analyze_overall_properties(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze the specific overall properties requested."""
        print("Analyzing overall properties...")

        stats = {}

        # _id analysis
        all_ids = [record.get("_id") for record in records if record.get("_id")]
        stats["_id"] = {
            "raw_count": len(all_ids),
            "unique_count": len(set(all_ids)),
            "duplicate_count": len(all_ids) - len(set(all_ids)),
        }

        # subject.id analysis
        subject_ids = [
            record.get("subject", {}).get("id")
            for record in records
            if record.get("subject", {}).get("id")
        ]
        stats["subject_id"] = {"raw_count": len(subject_ids), "unique_count": len(set(subject_ids))}

        # object.id analysis
        object_ids = [
            record.get("object", {}).get("id")
            for record in records
            if record.get("object", {}).get("id")
        ]
        stats["object_id"] = {"raw_count": len(object_ids), "unique_count": len(set(object_ids))}

        # subject.description analysis
        stats["subject_description"] = self._count_property_stats(records, "subject.description")

        # object.description analysis
        stats["object_description"] = self._count_property_stats(records, "object.description")

        # subject.xrefs analysis
        stats["subject_xrefs"] = self._count_property_stats(records, "subject.xrefs")

        # object.xrefs analysis
        stats["object_xrefs"] = self._count_property_stats(records, "object.xrefs")

        # association.publication analysis
        stats["association_publication"] = self._count_property_stats(
            records, "association.publication"
        )

        return stats

    def _analyze_relationship_type_properties(
            self, records: List[Dict[str, Any]], relationship_type: str
    ) -> Dict[str, Any]:
        """Analyze all properties for a specific relationship type."""
        print(f"Analyzing properties for {relationship_type}...")

        if not records:
            return {}

        subject_props = set()
        object_props = set()
        association_props = set()

        for record in records:
            if record.get("subject"):
                subject_props.update(self._extract_all_properties(record["subject"]))
            if record.get("object"):
                object_props.update(self._extract_all_properties(record["object"]))
            if record.get("association"):
                association_props.update(self._extract_all_properties(record["association"]))

        self.all_subject_properties.update(subject_props)
        self.all_object_properties.update(object_props)
        self.all_association_properties.update(association_props)

        stats = {
            "record_count": len(records),
            "subject_properties": {},
            "object_properties": {},
            "association_properties": {},
        }

        for prop in sorted(subject_props):
            stats["subject_properties"][prop] = self._count_property_stats(
                records, f"subject.{prop}"
            )

        for prop in sorted(object_props):
            stats["object_properties"][prop] = self._count_property_stats(records, f"object.{prop}")

        for prop in sorted(association_props):
            stats["association_properties"][prop] = self._count_property_stats(
                records, f"association.{prop}"
            )

        return stats

    def generate_comprehensive_stats(self, file_path: Optional[str] = None) -> Dict[str, Any]:
        """Generate comprehensive statistics for all records."""
        if file_path is None:
            file_path = self._find_latest_jsonl_file()
            if file_path is None:
                return {}

        print(f"Starting comprehensive analysis of: {file_path}")

        relationship_groups = defaultdict(list)
        all_records = []

        for record in self._load_jsonl_records(file_path):
            all_records.append(record)
            relationship_type = record.get("association", {}).get("category")
            if relationship_type:
                relationship_groups[relationship_type].append(record)

        if not all_records:
            print("No records found!")
            return {}

        print(f"Loaded {len(all_records)} total records")
        print(f"Found {len(relationship_groups)} relationship types:")
        for rel_type, records in relationship_groups.items():
            print(f"  - {rel_type}: {len(records)} records")

        overall_stats = self._analyze_overall_properties(all_records)

        relationship_stats = {}
        for rel_type, records in relationship_groups.items():
            relationship_stats[rel_type] = self._analyze_relationship_type_properties(
                records, rel_type
            )

        stats_report = {
            "metadata": {
                "analysis_date": datetime.now().isoformat(),
                "source_file": file_path,
                "total_records": len(all_records),
                "relationship_types": list(relationship_groups.keys()),
                "relationship_type_counts": {
                    rel_type: len(records) for rel_type, records in relationship_groups.items()
                },
            },
            "overall_statistics": overall_stats,
            "relationship_type_analysis": relationship_stats,
            "property_summary": {
                "all_subject_properties": sorted(list(self.all_subject_properties)),
                "all_object_properties": sorted(list(self.all_object_properties)),
                "all_association_properties": sorted(list(self.all_association_properties)),
            },
        }

        return stats_report

    def save_stats_report(
            self, stats: Dict[str, Any], output_filename: Optional[str] = None
    ) -> str:
        """Save the statistics report to JSON file."""
        if output_filename is None:
            timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
            output_filename = f"{timestamp}_hmdb_comprehensive_stats.json"

        report_path = os.path.join(self.report_dir, output_filename)

        with open(report_path, "w", encoding="utf-8") as f:
            json.dump(stats, f, indent=2, sort_keys=True, default=str)

        print(f"Statistics report saved to: {report_path}")
        return report_path

    def run_full_analysis(self, file_path: Optional[str] = None) -> Dict[str, Any]:
        """Run the complete analysis and save the report."""
        print("Starting full HMDB record analysis...")

        stats = self.generate_comprehensive_stats(file_path)

        if stats:
            report_path = self.save_stats_report(stats)
            print(f"Analysis complete! Report saved to: {report_path}")

            metadata = stats.get("metadata", {})
            print("\nSummary:")
            print(f"Total records analyzed: {metadata.get('total_records', 0)}")
            print(f"Relationship types found: {len(metadata.get('relationship_types', []))}")
            for rel_type, count in metadata.get("relationship_type_counts", {}).items():
                print(f"  - {rel_type}: {count}")
        else:
            print("Analysis failed - no data found.")

        return stats


if __name__ == "__main__":
    import sys

    file_path = sys.argv[1] if len(sys.argv) > 1 else None

    reporter = HMDBRecordStatsReporter()
    results = reporter.run_full_analysis(file_path)
