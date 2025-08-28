"""
Record deduplicator to remove duplicate records
from a list or iterator of records. It merges xrefs from duplicates
while preserving other fields from the first occurrence.
"""

import copy
import json
from typing import Dict, Iterator, List, Union


def _merge_xrefs(target_xrefs: Dict, source_xrefs: Dict) -> None:
    """
    Merges source_xrefs into target_xrefs in-place, combining values
    and ensuring no duplicates are created.
    Enhanced version with better type handling.
    """
    for key, source_value in source_xrefs.items():
        if key not in target_xrefs:
            target_xrefs[key] = source_value
            continue

        target_value = target_xrefs[key]

        combined_values = []
        if isinstance(target_value, list):
            combined_values.extend(target_value)
        else:
            combined_values.append(target_value)

        if isinstance(source_value, list):
            combined_values.extend(source_value)
        else:
            combined_values.append(source_value)

        unique_items = []
        seen_representations = set()
        for item in combined_values:
            if isinstance(item, dict):
                representation = json.dumps(item, sort_keys=True)
            elif isinstance(item, (list, tuple)):
                representation = json.dumps(item, sort_keys=True)
            else:
                representation = item

            if representation not in seen_representations:
                unique_items.append(item)
                seen_representations.add(representation)

        if len(unique_items) == 1:
            target_xrefs[key] = unique_items[0]
        else:
            target_xrefs[key] = unique_items


def _create_fingerprint(record: Dict) -> str:
    """
    Creates a unique, hashable fingerprint for a record by serializing its
    subject, object, and association nodes, excluding 'xrefs'.
    Enhanced to handle edge cases and provide more stable fingerprints.
    """
    fingerprint_data = {}

    for key in ["subject", "object", "association", "_id"]:
        if key in record:
            node_copy = copy.deepcopy(record[key])

            if key in ["subject", "object"] and isinstance(node_copy, dict):
                node_copy.pop("xrefs", None)

            fingerprint_data[key] = node_copy

    return json.dumps(fingerprint_data, sort_keys=True, default=str)


def _merge_records(target_record: Dict, source_record: Dict) -> None:
    """
    Merge source_record into target_record in-place.
    Focuses on merging xrefs while preserving other fields from target.
    """
    # merge subject xrefs
    if (
        "subject" in source_record
        and "xrefs" in source_record["subject"]
        and source_record["subject"]["xrefs"]
    ):

        if "subject" not in target_record:
            target_record["subject"] = {}
        if "xrefs" not in target_record["subject"]:
            target_record["subject"]["xrefs"] = {}

        _merge_xrefs(target_record["subject"]["xrefs"], source_record["subject"]["xrefs"])

    if (
        "object" in source_record
        and "xrefs" in source_record["object"]
        and source_record["object"]["xrefs"]
    ):

        if "object" not in target_record:
            target_record["object"] = {}
        if "xrefs" not in target_record["object"]:
            target_record["object"]["xrefs"] = {}

        _merge_xrefs(target_record["object"]["xrefs"], source_record["object"]["xrefs"])


def deduplicate_records_list(records: List[Dict]) -> List[Dict]:
    """
    Deduplicates a list of records (your original approach).
    Most memory-efficient for batch processing.
    """
    if not records:
        return []

    processed_records: Dict[str, Dict] = {}

    for record in records:
        fingerprint = _create_fingerprint(record)

        if fingerprint not in processed_records:
            processed_records[fingerprint] = json.loads(json.dumps(record))
        else:
            _merge_records(processed_records[fingerprint], record)

    return list(processed_records.values())


def deduplicate_records_iterator(records: Iterator[Dict]) -> Iterator[Dict]:
    """
    Deduplicates records from an iterator (streaming approach).
    Memory usage grows with unique records, but processes one at a time.
    """
    processed_records: Dict[str, Dict] = {}

    for record in records:
        fingerprint = _create_fingerprint(record)

        if fingerprint not in processed_records:
            processed_records[fingerprint] = json.loads(json.dumps(record))
        else:
            _merge_records(processed_records[fingerprint], record)

    for record in processed_records.values():
        yield record


def deduplicate_records(
    records: Union[List[Dict], Iterator[Dict]]
) -> Union[List[Dict], Iterator[Dict]]:
    """
    Smart deduplication that adapts to input type.

    Args:
        records: Either a list or iterator of record dictionaries

    Returns:
        Deduplicated records in the same format as input

    Examples:
        # List input -> List output
        clean_records = deduplicate_records(record_list)

        # Iterator input -> Iterator output
        clean_records = deduplicate_records(parser.parse_associations())
    """
    if isinstance(records, list):
        return deduplicate_records_list(records)
    else:
        return deduplicate_records_iterator(records)


def deduplicate_with_stats(records: Union[List[Dict], Iterator[Dict]]) -> Dict:
    """
    Deduplicates records and returns detailed statistics about the process.

    Returns:
        {
            'deduplicated_records': List[Dict],
            'stats': {
                'original_count': int,
                'deduplicated_count': int,
                'duplicates_removed': int,
                'records_with_merged_xrefs': int,
                'duplicate_rate': float
            }
        }
    """
    if not isinstance(records, list):
        records = list(records)

    original_count = len(records)
    processed_records: Dict[str, Dict] = {}
    records_with_merged_xrefs = 0

    for record in records:
        fingerprint = _create_fingerprint(record)

        if fingerprint not in processed_records:
            processed_records[fingerprint] = json.loads(json.dumps(record))
        else:
            existing_record = processed_records[fingerprint]

            will_merge_xrefs = False
            for node_key in ["subject", "object"]:
                if node_key in record and "xrefs" in record[node_key] and record[node_key]["xrefs"]:
                    will_merge_xrefs = True
                    break

            if will_merge_xrefs:
                records_with_merged_xrefs += 1

            _merge_records(existing_record, record)

    deduplicated_records = list(processed_records.values())
    deduplicated_count = len(deduplicated_records)

    return {
        "deduplicated_records": deduplicated_records,
        "stats": {
            "original_count": original_count,
            "deduplicated_count": deduplicated_count,
            "duplicates_removed": original_count - deduplicated_count,
            "records_with_merged_xrefs": records_with_merged_xrefs,
            "duplicate_rate": (original_count - deduplicated_count) / original_count * 100
            if original_count > 0
            else 0.0,
        },
    }


def deduplicate_and_merge(records: List[Dict]) -> List[Dict]:
    """Backwards compatible alias for your existing code."""
    return deduplicate_records_list(records)


class DeduplicationIntegrationExamples:
    """Examples of how to integrate into existing systems."""

    @staticmethod
    def record_manager_batch_approach(metabolite_parser):
        """Batch processing - good for smaller datasets."""
        raw_records = list(metabolite_parser.parse_microbe_metabolite())
        return deduplicate_records(raw_records)

    @staticmethod
    def record_manager_streaming_approach(metabolite_parser):
        """Streaming approach - memory efficient for large datasets."""
        return list(deduplicate_records(metabolite_parser.parse_microbe_metabolite()))

    @staticmethod
    def record_manager_with_stats(metabolite_parser):
        """Get deduplication statistics for monitoring."""
        raw_records = list(metabolite_parser.parse_microbe_metabolite())
        result = deduplicate_with_stats(raw_records)

        print(
            f"Deduplicated {result['stats']['duplicates_removed']} duplicates "
            f"({result['stats']['duplicate_rate']:.1f}% duplicate rate)"
        )

        return result["deduplicated_records"]


def test_hybrid_deduplication():
    """Test the hybrid deduplication with various scenarios."""

    test_records = [
        # Exact duplicate
        {
            "_id": "test1",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"uniprot": "P1"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"kegg": "K1"}},
            "association": {"type": "participates_in"},
        },
        # Exact duplicate of above
        {
            "_id": "test1",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"uniprot": "P1"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"kegg": "K1"}},
            "association": {"type": "participates_in"},
        },
        # Same relationship, additional xrefs
        {
            "_id": "test1",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"uniprot": "P1", "pdb": "1ABC"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"kegg": "K1", "smpdb": "S1"}},
            "association": {"type": "participates_in"},
        },
        # Completely different record
        {
            "_id": "test2",
            "subject": {"id": "C", "name": "protein_c"},
            "object": {"id": "D", "name": "pathway_d"},
            "association": {"type": "involved_in"},
        },
    ]

    # Test batch deduplication
    result = deduplicate_with_stats(test_records)

    assert len(result["deduplicated_records"]) == 2, "Should have 2 unique records"
    assert result["stats"]["duplicates_removed"] == 2, "Should remove 2 duplicates"

    # Check that xrefs were properly merged
    merged_record = next(r for r in result["deduplicated_records"] if r["_id"] == "test1")
    assert "pdb" in merged_record["subject"]["xrefs"], "Should have merged PDB xref"
    assert "smpdb" in merged_record["object"]["xrefs"], "Should have merged SMPDB xref"

    print("Deduplication test passed!")
    return True


if __name__ == "__main__":
    test_hybrid_deduplication()
