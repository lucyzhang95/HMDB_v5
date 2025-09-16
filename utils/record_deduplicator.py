"""
Two-stage record deduplicator:
1. Merge xrefs when non-xref fields are identical
2. Remove exact duplicates after merging
"""

import copy
import json
from typing import Dict, Iterator, List, Union


def _merge_xrefs(target_xrefs: Dict, source_xrefs: Dict) -> None:
    """
    Merges source_xrefs into target_xrefs in-place, combining values
    and ensuring no duplicates are created.
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


def _create_merge_fingerprint(record: Dict) -> str:
    """
    Creates fingerprint for merging records with same non-xref fields.
    Excludes xrefs from subject/object and pmid from association.publication.
    """
    fingerprint_data = {}

    for key in ["subject", "object", "association", "_id"]:
        if key in record:
            node_copy = copy.deepcopy(record[key])

            # Remove xrefs for fingerprint calculation
            if key in ["subject", "object"] and isinstance(node_copy, dict):
                node_copy.pop("xrefs", None)

            # Remove pmid from association.publication for fingerprint calculation
            elif key == "association" and isinstance(node_copy, dict):
                if "publication" in node_copy and isinstance(node_copy["publication"], dict):
                    publication_copy = node_copy["publication"].copy()
                    publication_copy.pop("pmid", None)
                    node_copy["publication"] = publication_copy

            fingerprint_data[key] = node_copy

    return json.dumps(fingerprint_data, sort_keys=True, default=str)


def _create_exact_fingerprint(record: Dict) -> str:
    """
    Creates fingerprint for exact duplicate detection.
    Includes ALL fields including merged xrefs.
    """
    return json.dumps(record, sort_keys=True, default=str)


def _merge_pmids(target_pmids: List, source_pmids: List) -> List:
    """Merge two lists of PMIDs, removing duplicates while preserving order."""
    combined_pmids = target_pmids + source_pmids
    unique_pmids = []
    seen = set()

    for pmid in combined_pmids:
        if pmid not in seen:
            unique_pmids.append(pmid)
            seen.add(pmid)

    return unique_pmids


def _has_better_evidence(record: Dict) -> bool:
    """
    Check if a record has better evidence quality.
    Returns True if record has publication and non-zero evidence code.
    """
    has_publication = (
            "association" in record
            and "publication" in record["association"]
            and record["association"]["publication"] is not None
    )
    has_good_evidence = (
            "association" in record
            and "has_evidence" in record["association"]
            and record["association"]["has_evidence"] != "ECO:0000000"
    )

    return has_publication and has_good_evidence


def _merge_records(target_record: Dict, source_record: Dict) -> None:
    """
    Merge source_record into target_record in-place.
    Focuses on merging xrefs and pmids while preserving other fields from target.
    If source has better evidence (publication + non-zero ECO), replace target entirely.
    """
    if _has_better_evidence(source_record) and not _has_better_evidence(target_record):
        # replace target with source (source has better evidence)
        target_record.clear()
        target_record.update(copy.deepcopy(source_record))
        return
    elif _has_better_evidence(target_record) and not _has_better_evidence(source_record):
        return

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

    # merge object xrefs
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

    # merge association publication pmids
    if (
            "association" in source_record
            and "publication" in source_record["association"]
            and isinstance(source_record["association"]["publication"], dict)
            and "pmid" in source_record["association"]["publication"]
    ):
        if "association" not in target_record:
            target_record["association"] = {}
        if "publication" not in target_record["association"]:
            target_record["association"]["publication"] = {}

        source_pub = source_record["association"]["publication"]
        target_pub = target_record["association"]["publication"]

        # check if publications match (excluding pmid field)
        source_pub_key = {k: v for k, v in source_pub.items() if k != "pmid"}
        target_pub_key = {k: v for k, v in target_pub.items() if k != "pmid"}

        if source_pub_key == target_pub_key or not target_pub_key:
            source_pmids = source_pub["pmid"]
            target_pmids = target_pub.get("pmid", [])

            if isinstance(source_pmids, str):
                source_pmids = [source_pmids]
            if isinstance(target_pmids, str):
                target_pmids = [target_pmids]
            elif not isinstance(target_pmids, list):
                target_pmids = []

            # merge pmids
            merged_pmids = _merge_pmids(target_pmids, source_pmids)

            target_pub.update(source_pub_key)
            target_pub["pmid"] = merged_pmids


def deduplicate_records(records: Union[List[Dict], Iterator[Dict]]) -> List[Dict]:
    """
    Two-stage deduplication:
    1. Merge xrefs for records with identical non-xref fields
    2. Remove exact duplicates after merging

    Args:
        records: List or iterator of record dictionaries

    Returns:
        List of deduplicated records
    """
    if not isinstance(records, list):
        records = list(records)

    if not records:
        return []

    # merge records with same non-xref fields
    merge_groups = {}

    for record in records:
        merge_fingerprint = _create_merge_fingerprint(record)

        if merge_fingerprint not in merge_groups:
            merge_groups[merge_fingerprint] = json.loads(json.dumps(record))
        else:
            _merge_records(merge_groups[merge_fingerprint], record)

    merged_records = list(merge_groups.values())

    # remove exact duplicates from merged records
    exact_fingerprints = set()
    final_records = []

    for record in merged_records:
        exact_fingerprint = _create_exact_fingerprint(record)

        if exact_fingerprint not in exact_fingerprints:
            exact_fingerprints.add(exact_fingerprint)
            final_records.append(record)

    return final_records


def deduplicate_with_stats(records: Union[List[Dict], Iterator[Dict]]) -> Dict:
    """
    Two-stage deduplication with statistics.

    Returns:
        {
            'deduplicated_records': List[Dict],
            'stats': {
                'original_count': int,
                'after_merge_count': int,
                'final_count': int,
                'records_merged': int,
                'exact_duplicates_removed': int,
                'duplicate_rate': float
            }
        }
    """
    if not isinstance(records, list):
        records = list(records)

    original_count = len(records)

    if not records:
        return {
            "deduplicated_records": [],
            "stats": {
                "original_count": 0,
                "after_merge_count": 0,
                "final_count": 0,
                "records_merged": 0,
                "exact_duplicates_removed": 0,
                "duplicate_rate": 0.0
            }
        }

    # merge records with same non-xref fields
    merge_groups = {}
    records_merged = 0

    for record in records:
        merge_fingerprint = _create_merge_fingerprint(record)

        if merge_fingerprint not in merge_groups:
            merge_groups[merge_fingerprint] = json.loads(json.dumps(record))
        else:
            _merge_records(merge_groups[merge_fingerprint], record)
            records_merged += 1

    merged_records = list(merge_groups.values())
    after_merge_count = len(merged_records)

    # remove exact duplicates from merged records
    exact_fingerprints = set()
    final_records = []
    exact_duplicates_removed = 0

    for record in merged_records:
        exact_fingerprint = _create_exact_fingerprint(record)

        if exact_fingerprint not in exact_fingerprints:
            exact_fingerprints.add(exact_fingerprint)
            final_records.append(record)
        else:
            exact_duplicates_removed += 1

    final_count = len(final_records)

    stats = {
        "original_count": original_count,
        "after_merge_count": after_merge_count,
        "final_count": final_count,
        "records_merged": records_merged,
        "exact_duplicates_removed": exact_duplicates_removed,
        "duplicate_rate": (original_count - final_count) / original_count * 100 if original_count > 0 else 0.0
    }

    return {
        "deduplicated_records": final_records,
        "stats": stats
    }


def test_two_stage_deduplication():
    """Test the deduplication process."""

    test_records = [
        # Records 1-2: same non-xref fields, different subject and object xrefs as well as different pmid
        {
            "_id": "A_participates_in_B",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"uniprot": "P1"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"kegg": "K1"}},
            "association": {
                "category": "participates_in",
                "publication": {
                    "category": "biolink:Publication",
                    "pmid": ["PMID:12345", "PMID:67890"],
                },
            },
        },
        {
            "_id": "A_participates_in_B",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"pdb": "1ABC"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"smpdb": "S1"}},
            "association": {
                "category": "participates_in",
                "publication": {"category": "biolink:Publication", "pmid": "PMID:11111"},
            },
        },
        # Record 3: Exact duplicate after merging -> should be removed in stage 2
        {
            "_id": "A_participates_in_B",
            "subject": {"id": "A", "name": "protein_a", "xrefs": {"uniprot": "P1", "pdb": "1ABC"}},
            "object": {"id": "B", "name": "pathway_b", "xrefs": {"kegg": "K1", "smpdb": "S1"}},
            "association": {
                "category": "participates_in",
                "publication": {
                    "category": "biolink:Publication",
                    "pmid": ["PMID:12345", "PMID:67890", "PMID:11111"],
                },
            },
        },
        # Record 4: completely different -> should remain separate
        {
            "_id": "C_involved_in_D",
            "subject": {"id": "C", "name": "protein_c"},
            "object": {"id": "D", "name": "pathway_d"},
            "association": {"category": "involved_in"},
        },
    ]

    print("=== TWO-STAGE DEDUPLICATION TEST ===")
    print(f"Original records: {len(test_records)}")

    result = deduplicate_with_stats(test_records)
    stats = result["stats"]
    records = result["deduplicated_records"]

    print(f"After merge stage: {stats['after_merge_count']}")
    print(f"Records merged in stage 1: {stats['records_merged']}")
    print(f"Exact duplicates removed in stage 2: {stats['exact_duplicates_removed']}")
    print(f"Final unique records: {stats['final_count']}")
    print(f"Overall duplicate rate: {stats['duplicate_rate']:.1f}%")

    # Verify the merged record has all xrefs and pmids
    test1_record = next(r for r in records if r["_id"] == "A_participates_in_B")
    subject_xrefs = test1_record["subject"]["xrefs"]
    object_xrefs = test1_record["object"]["xrefs"]

    print(f"\nMerged subject xrefs: {list(subject_xrefs.keys())}")
    print(f"Merged object xrefs: {list(object_xrefs.keys())}")

    # check merged pmids
    if "association" in test1_record and "publication" in test1_record["association"]:
        pub = test1_record["association"]["publication"]
        if "pmid" in pub:
            print(f"Merged pmids: {pub['pmid']}")

    expected_subject_xrefs = {"uniprot", "pdb"}
    expected_object_xrefs = {"kegg", "smpdb"}
    expected_pmids = ["PMID:12345", "PMID:67890", "PMID:11111"]

    assert set(subject_xrefs.keys()) == expected_subject_xrefs
    assert set(object_xrefs.keys()) == expected_object_xrefs
    assert len(records) == 2

    # Check pmids were merged correctly
    if "association" in test1_record and "publication" in test1_record["association"]:
        pub = test1_record["association"]["publication"]
        assert "pmid" in pub, "Should have pmid field"
        assert pub["pmid"] == expected_pmids, f"Expected {expected_pmids}, got {pub['pmid']}"

    print("\n✅ Two-stage deduplication with PMID merging test passed!")
    return records


if __name__ == "__main__":
    test_two_stage_deduplication()
