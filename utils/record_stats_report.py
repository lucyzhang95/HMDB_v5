import json
import os
import pickle
import random
from collections import Counter, defaultdict
from datetime import datetime
from statistics import mean, median
from typing import Any, Dict, List


class HMDBRecordStatsReporter:
    """Generates comprehensive statistics for HMDB parsed records."""

    def __init__(self, report_dir="../reports"):
        self.report_dir = report_dir
        os.makedirs(self.report_dir, exist_ok=True)
        print("▶️ HMDBRecordStatsReporter initialized.")

    def _safe_get_numeric_stats(self, values: List[Any]) -> Dict[str, Any]:
        """Safely calculate numeric statistics for a list of values."""
        numeric_values = []
        for v in values:
            if v is not None:
                try:
                    if isinstance(v, dict):
                        if "average_molecular_weight" in v:
                            numeric_values.append(float(v["average_molecular_weight"]))
                        elif "monoisotopic_molecular_weight" in v:
                            numeric_values.append(float(v["monoisotopic_molecular_weight"]))
                    else:
                        numeric_values.append(float(v))
                except (ValueError, TypeError):
                    continue

        if not numeric_values:
            return {"count": 0, "min": None, "max": None, "mean": None, "median": None}

        return {
            "count": len(numeric_values),
            "min": min(numeric_values),
            "max": max(numeric_values),
            "mean": round(mean(numeric_values), 3),
            "median": round(median(numeric_values), 3),
        }

    def _extract_curie_prefix(self, curie: str) -> str:
        """Extract prefix from CURIE format ID."""
        if ":" in curie:
            return curie.split(":", 1)[0]
        return "unknown"

    def _count_xrefs(self, node: Dict[str, Any]) -> Dict[str, int]:
        """Count xref types in a node."""
        xrefs = node.get("xrefs", {})
        if not xrefs:
            return {}

        xref_counts = {}
        for key, value in xrefs.items():
            if value:  # only non-empty xrefs
                if isinstance(value, list):
                    xref_counts[key] = len(value)
                else:
                    xref_counts[key] = 1
        return xref_counts

    def _collect_unique_xrefs(self, node: Dict[str, Any]) -> Dict[str, set]:
        """Collect unique xref values by type from a node."""
        xrefs = node.get("xrefs", {})
        if not xrefs:
            return {}

        unique_xrefs = {}
        for key, value in xrefs.items():
            if value:  # only non-empty xrefs
                if isinstance(value, list):
                    processed_values = []
                    for item in value:
                        if isinstance(item, dict):
                            # Extract id or name from dict
                            if "id" in item:
                                processed_values.append(item["id"])
                            elif "name" in item:
                                processed_values.append(item["name"])
                            else:
                                processed_values.append(str(item))
                        else:
                            processed_values.append(str(item))
                    unique_xrefs[key] = set(processed_values)
                elif isinstance(value, dict):
                    if "id" in value:
                        unique_xrefs[key] = {value["id"]}
                    elif "name" in value:
                        unique_xrefs[key] = {value["name"]}
                    else:
                        unique_xrefs[key] = {str(value)}
                else:
                    unique_xrefs[key] = {str(value)}
        return unique_xrefs

    def _analyze_molecular_weight(self, nodes: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze molecular weight data from nodes."""
        mw_data = {"average": [], "monoisotopic": []}

        for node in nodes:
            mw = node.get("molecular_weight", {})
            if isinstance(mw, dict):
                if "average_molecular_weight" in mw and mw["average_molecular_weight"] is not None:
                    mw_data["average"].append(mw["average_molecular_weight"])
                if (
                    "monoisotopic_molecular_weight" in mw
                    and mw["monoisotopic_molecular_weight"] is not None
                ):
                    mw_data["monoisotopic"].append(mw["monoisotopic_molecular_weight"])

        return {
            "average": self._safe_get_numeric_stats(mw_data["average"]),
            "monoisotopic": self._safe_get_numeric_stats(mw_data["monoisotopic"]),
        }

    def _count_publication_pmids(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Count publication records and unique PMIDs."""
        records_with_pub = 0
        all_pmids = []

        for record in records:
            pub = record.get("association", {}).get("publication", {})
            if pub and pub.get("pmid"):
                records_with_pub += 1
                pmid = pub["pmid"]
                if isinstance(pmid, list):
                    all_pmids.extend(pmid)
                else:
                    all_pmids.append(pmid)

        unique_pmids = len(set(all_pmids))

        return {
            "records_with_publication": records_with_pub,
            "unique_pmid_count": unique_pmids,
            "total_pmid_count": len(all_pmids),
        }

    def _count_unique_items_in_list_field(
        self, records: List[Dict[str, Any]], field_path: List[str]
    ) -> Dict[str, int]:
        """Count records with field and unique items in list fields."""
        records_with_field = 0
        all_items = []

        for record in records:
            current = record
            for path_part in field_path:
                current = current.get(path_part, {})
                if not current:
                    break

            if current and isinstance(current, list) and len(current) > 0:
                records_with_field += 1
                for item in current:
                    if isinstance(item, dict):
                        if "id" in item:
                            all_items.append(item["id"])
                        elif "original_name" in item:
                            all_items.append(item["original_name"])
                        else:
                            all_items.append(str(item))
                    else:
                        all_items.append(item)

        return {"record_count": records_with_field, "unique_count": len(set(all_items))}

    def export_duplicate_records(self) -> Dict[str, Any]:
        """Export complete duplicate records to JSON file for investigation."""
        print("\n📍 Identifying and exporting duplicate records...")

        filepath = os.path.join("..", "cache", "hmdb_v5_parsed_records.pkl")

        try:
            with open(filepath, "rb") as in_f:
                combined_data = pickle.load(in_f)
                if not combined_data:
                    print("❌ No data loaded from pickle file.")
                    return {}
        except Exception as e:
            print(f"‼️ Error loading hmdb_v5_parsed_records.pkl: {e}")
            return {}

        relationship_types = [
            "microbe-metabolite",
            "metabolite-disease",
            "metabolite-protein",
            "metabolite-pathway",
            "protein-pathway",
            "protein-biological_process",
        ]

        all_records_with_type = []
        for rel_type in relationship_types:
            if rel_type in combined_data:
                for record in combined_data[rel_type]:
                    record_copy = record.copy()
                    record_copy["_relationship_type"] = rel_type
                    all_records_with_type.append(record_copy)

        print(f"Total records collected: {len(all_records_with_type)}")

        # find duplicates by _id
        id_to_records = defaultdict(list)
        for record in all_records_with_type:
            record_id = record.get("_id")
            if record_id:
                id_to_records[record_id].append(record)

        duplicate_records = {}
        duplicate_summary = {}

        for record_id, records in id_to_records.items():
            if len(records) > 1:
                duplicate_records[record_id] = {
                    "count": len(records),
                    "relationship_types": list(
                        set(r.get("_relationship_type", "unknown") for r in records)
                    ),
                    "records": records,
                }
                duplicate_summary[record_id] = {
                    "count": len(records),
                    "relationship_types": duplicate_records[record_id]["relationship_types"],
                }

        print(f"✅ Found {len(duplicate_records)} duplicate IDs")

        duplicate_stats = self._generate_duplicate_statistics(duplicate_records)

        export_data = {
            "metadata": {
                "export_date": datetime.now().isoformat(),
                "total_duplicate_ids": len(duplicate_records),
                "total_duplicate_records": sum(
                    data["count"] for data in duplicate_records.values()
                ),
                "source_file": filepath,
                "relationship_types_analyzed": relationship_types,
            },
            "duplicate_summary": duplicate_summary,
            "duplicate_records": duplicate_records,
            "statistics": duplicate_stats,
        }

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        export_path = os.path.join(self.report_dir, f"duplicate_records_{timestamp}.json")

        with open(export_path, "w") as f:
            json.dump(export_data, f, indent=2, sort_keys=True, default=str)

        print(f"✅ Duplicate records exported to: {export_path}")
        print(f"✅ Exported {len(duplicate_records)} duplicate record groups")

        return export_data

    def _generate_duplicate_statistics(self, duplicate_records: Dict[str, Any]) -> Dict[str, Any]:
        """Generate statistics about the duplicate records."""
        if not duplicate_records:
            return {}

        count_distribution = Counter()
        relationship_type_pairs = Counter()
        cross_relationship_duplicates = 0

        for _, data in duplicate_records.items():
            count = data["count"]
            rel_types = data["relationship_types"]

            count_distribution[count] += 1

            if len(rel_types) > 1:
                cross_relationship_duplicates += 1
                rel_types_sorted = tuple(sorted(rel_types))
                relationship_type_pairs[rel_types_sorted] += 1

        most_common_counts = dict(count_distribution.most_common(10))
        most_common_rel_pairs = dict(relationship_type_pairs.most_common(10))

        return {
            "count_distribution": most_common_counts,
            "cross_relationship_duplicates": cross_relationship_duplicates,
            "relationship_type_pairs": most_common_rel_pairs,
            "max_duplicates_for_single_id": max(count_distribution.keys())
            if count_distribution
            else 0,
            "total_duplicate_groups": len(duplicate_records),
        }

    def generate_record_stats(self) -> Dict[str, Any]:
        """Generate comprehensive statistics for HMDB records."""
        print("\n⚙️ Generating HMDB record statistics...")
        filepath = os.path.join("..", "cache", "hmdb_v5_parsed_records.pkl")

        try:
            with open(filepath, "rb") as in_f:
                combined_data = pickle.load(in_f)
                if combined_data:
                    print("✅ Loaded data from hmdb_v5_parsed_records.pkl")
                else:
                    print("❌ File exists but no data loaded.")
                    return {}
        except Exception as e:
            print(f"‼️ Error loading hmdb_v5_parsed_records.pkl: {e}")
            return {}

        # HMDB 6 relationship types
        relationship_types = [
            "microbe-metabolite",
            "metabolite-disease",
            "metabolite-protein",
            "metabolite-pathway",
            "protein-pathway",
            "protein-biological_process",
        ]

        stats = {
            "data_source": "HMDB",
            "total_relationships": len(
                [k for k in combined_data.keys() if k in relationship_types]
            ),
            "relationship_types": [k for k in relationship_types if k in combined_data],
            "metadata": {
                "report_version": "1.0",
                "analysis_date": datetime.now().isoformat(),
            },
        }

        # total record counts
        total_records = 0
        relationship_counts = {}
        for rel_type in relationship_types:
            if rel_type in combined_data:
                count = len(combined_data[rel_type])
                relationship_counts[rel_type] = count
                total_records += count

        stats["metadata"]["total_records_analyzed"] = total_records
        stats["metadata"]["total_records_analyzed_by_relationship"] = relationship_counts

        # each relationship type
        relationship_stats = {}

        for rel_type in relationship_types:
            if rel_type in combined_data:
                records = combined_data[rel_type]
                rel_stats = self._analyze_relationship(rel_type, records)
                relationship_stats[rel_type] = rel_stats

        stats["relationship_analysis"] = relationship_stats

        stats.update(self._calculate_overall_stats(combined_data, relationship_types))

        print(f"Generated comprehensive statistics for {total_records} total records.")
        return stats

    def _calculate_overall_stats(
        self, combined_data: Dict, relationship_types: List[str]
    ) -> Dict[str, Any]:
        """Calculate overall statistics across all relationships."""
        all_records = []
        for rel_type in relationship_types:
            if rel_type in combined_data:
                all_records.extend(combined_data[rel_type])

        # all evidence types
        evidence_counter = Counter()
        for record in all_records:
            evidence_type = record.get("association", {}).get("evidence_type")
            if evidence_type:
                if isinstance(evidence_type, list):
                    evidence_counter.update(evidence_type)
                else:
                    evidence_counter[evidence_type] += 1

        # all ID duplication check
        all_record_ids = [record.get("_id") for record in all_records]
        id_counter = Counter(all_record_ids)
        duplicates = {id_val: count for id_val, count in id_counter.items() if count > 1}

        count_groups = defaultdict(list)
        for id_val, count in id_counter.items():
            if count > 1:
                count_groups[count].append(id_val)

        sampled_ids_by_count = {}
        for count, id_list in count_groups.items():
            sample_size = min(3, len(id_list))
            sampled_ids_by_count[count] = random.sample(id_list, sample_size)

        # all CURIE analysis
        overall_subject_curies = []
        overall_object_curies = []
        for record in all_records:
            subject_id = record.get("subject", {}).get("id", "")
            object_id = record.get("object", {}).get("id", "")
            if subject_id:
                overall_subject_curies.append(self._extract_curie_prefix(subject_id))
            if object_id:
                overall_object_curies.append(self._extract_curie_prefix(object_id))

        # all description analysis
        subject_desc_count = sum(
            1 for record in all_records if record.get("subject", {}).get("description")
        )
        object_desc_count = sum(
            1 for record in all_records if record.get("object", {}).get("description")
        )

        # all xrefs analysis
        overall_subject_xrefs = defaultdict(int)
        overall_object_xrefs = defaultdict(int)
        overall_subject_unique_xrefs = defaultdict(set)
        overall_object_unique_xrefs = defaultdict(set)

        for record in all_records:
            subject_xrefs = self._count_xrefs(record.get("subject", {}))
            object_xrefs = self._count_xrefs(record.get("object", {}))
            subject_unique_xrefs = self._collect_unique_xrefs(record.get("subject", {}))
            object_unique_xrefs = self._collect_unique_xrefs(record.get("object", {}))

            for xref_type, count in subject_xrefs.items():
                overall_subject_xrefs[xref_type] += count
            for xref_type, count in object_xrefs.items():
                overall_object_xrefs[xref_type] += count

            for xref_type, unique_values in subject_unique_xrefs.items():
                overall_subject_unique_xrefs[xref_type].update(unique_values)
            for xref_type, unique_values in object_unique_xrefs.items():
                overall_object_unique_xrefs[xref_type].update(unique_values)

        # all publication analysis
        overall_pub_stats = self._count_publication_pmids(all_records)

        return {
            "overall_evidence_types": dict(evidence_counter),
            "overall_id_analysis": {
                "total_ids": len(all_record_ids),
                "unique_ids": len(id_counter),
                "duplicate_count": len(duplicates),
                "duplicates": sampled_ids_by_count,
                "duplicates_distribution": dict(Counter(k for k in id_counter.values() if k > 1)),
            },
            "overall_subject_curie_stats": dict(Counter(overall_subject_curies)),
            "overall_unique_subject_curie_count": len(set(overall_subject_curies)),
            "overall_object_curie_stats": dict(Counter(overall_object_curies)),
            "overall_unique_object_curie_count": len(set(overall_object_curies)),
            "overall_description_stats": {
                "subject_descriptions": subject_desc_count,
                "subject_description_percentage": round(
                    (subject_desc_count / len(all_records)) * 100, 2
                ),
                "object_descriptions": object_desc_count,
                "object_description_percentage": round(
                    (object_desc_count / len(all_records)) * 100, 2
                ),
            },
            "overall_xref_stats": {
                "subject_xrefs": dict(overall_subject_xrefs),
                "object_xrefs": dict(overall_object_xrefs),
                "subject_unique_xrefs": {
                    k: len(v) for k, v in overall_subject_unique_xrefs.items()
                },
                "object_unique_xrefs": {k: len(v) for k, v in overall_object_unique_xrefs.items()},
            },
            "overall_publication_stats": overall_pub_stats,
        }

    def _analyze_relationship(self, rel_type: str, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze statistics for a specific relationship type."""
        print(f"-> Analyzing {rel_type} relationship ({len(records)} records)...")

        rel_stats = {
            "record_count": len(records),
            "evidence_types": {},
            "subject_curie_stats": {},
            "object_curie_stats": {},
            "subject_description_count": 0,
            "object_description_count": 0,
            "subject_xref_stats": defaultdict(int),
            "object_xref_stats": defaultdict(int),
            "subject_unique_xref_stats": defaultdict(set),
            "object_unique_xref_stats": defaultdict(set),
            "publication_stats": {},
        }

        evidence_counter = Counter()
        subject_curie_counter = Counter()
        object_curie_counter = Counter()

        for record in records:
            # evidence types
            evidence_type = record.get("association", {}).get("evidence_type")
            if evidence_type:
                if isinstance(evidence_type, list):
                    evidence_counter.update(evidence_type)
                else:
                    evidence_counter[evidence_type] += 1

            # CURIE prefixes
            subject_id = record.get("subject", {}).get("id", "")
            object_id = record.get("object", {}).get("id", "")

            if subject_id:
                subject_curie_counter[self._extract_curie_prefix(subject_id)] += 1
            if object_id:
                object_curie_counter[self._extract_curie_prefix(object_id)] += 1

            # descriptions
            if record.get("subject", {}).get("description"):
                rel_stats["subject_description_count"] += 1
            if record.get("object", {}).get("description"):
                rel_stats["object_description_count"] += 1

            # xrefs
            subject_xrefs = self._count_xrefs(record.get("subject", {}))
            object_xrefs = self._count_xrefs(record.get("object", {}))
            subject_unique_xrefs = self._collect_unique_xrefs(record.get("subject", {}))
            object_unique_xrefs = self._collect_unique_xrefs(record.get("object", {}))

            for xref_type, count in subject_xrefs.items():
                rel_stats["subject_xref_stats"][xref_type] += count
            for xref_type, count in object_xrefs.items():
                rel_stats["object_xref_stats"][xref_type] += count

            for xref_type, unique_values in subject_unique_xrefs.items():
                rel_stats["subject_unique_xref_stats"][xref_type].update(unique_values)
            for xref_type, unique_values in object_unique_xrefs.items():
                rel_stats["object_unique_xref_stats"][xref_type].update(unique_values)

        rel_stats["evidence_types"] = dict(evidence_counter)
        rel_stats["subject_curie_stats"] = dict(subject_curie_counter)
        rel_stats["object_curie_stats"] = dict(object_curie_counter)
        rel_stats["subject_xref_stats"] = dict(rel_stats["subject_xref_stats"])
        rel_stats["object_xref_stats"] = dict(rel_stats["object_xref_stats"])
        rel_stats["subject_unique_xref_stats"] = {
            k: len(v) for k, v in rel_stats["subject_unique_xref_stats"].items()
        }
        rel_stats["object_unique_xref_stats"] = {
            k: len(v) for k, v in rel_stats["object_unique_xref_stats"].items()
        }

        # publication stats
        rel_stats["publication_stats"] = self._count_publication_pmids(records)

        # relationship-specific analysis
        if rel_type == "microbe-metabolite":
            rel_stats.update(self._analyze_microbe_metabolite_specific(records))
        elif rel_type == "metabolite-disease":
            rel_stats.update(self._analyze_metabolite_disease_specific(records))
        elif rel_type == "metabolite-protein":
            rel_stats.update(self._analyze_metabolite_protein_specific(records))
        elif rel_type == "protein-pathway":
            rel_stats.update(self._analyze_protein_pathway_specific(records))
        elif rel_type == "protein-biological_process":
            rel_stats.update(self._analyze_protein_biological_process_specific(records))

        return rel_stats

    def _analyze_microbe_metabolite_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze microbe-metabolite specific statistics."""
        organism_types = []
        ranks = []
        object_nodes = []
        logp_values = []
        melting_point_count = 0

        for record in records:
            subject = record.get("subject", {})
            object_node = record.get("object", {})
            object_nodes.append(object_node)

            # organism type and rank
            organism_type = subject.get("organism_type")
            if organism_type:
                organism_types.append(organism_type)

            rank = subject.get("rank")
            if rank:
                ranks.append(rank)

            # object properties
            logp = object_node.get("logp")
            if logp is not None:
                logp_values.append(logp)

            if object_node.get("melting_point"):
                melting_point_count += 1

        return {
            "organism_types": dict(Counter(organism_types)),
            "ranks": dict(Counter(ranks)),
            "object_molecular_weight_stats": self._analyze_molecular_weight(object_nodes),
            "object_logp_stats": self._safe_get_numeric_stats(logp_values),
            "melting_point_count": melting_point_count,
            "cellular_component_stats": self._count_unique_items_in_list_field(
                records, ["object", "cellular_component"]
            ),
            "biosample_stats": self._count_unique_items_in_list_field(
                records, ["object", "biosample"]
            ),
            "anatomical_entity_stats": self._count_unique_items_in_list_field(
                records, ["object", "anatomical_entity"]
            ),
        }

    def _analyze_metabolite_disease_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze metabolite-disease specific statistics."""
        return {}

    def _analyze_metabolite_protein_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze metabolite-protein specific statistics."""
        protein_types = []

        for record in records:
            object_node = record.get("object", {})
            protein_type = object_node.get("protein_type")
            if protein_type:
                protein_types.append(protein_type)

        return {
            "metabolite_protein_object_protein_type_stats": {
                "record_count": len(protein_types),
                "unique_count": len(set(protein_types)),
                "distribution": dict(Counter(protein_types)),
            }
        }

    def _analyze_protein_pathway_specific(self, records: List[Dict[str, Any]]) -> Dict[str, Any]:
        """Analyze protein-pathway specific statistics."""
        functions = []
        molecular_weights = []
        transmembrane_count = 0
        signal_region_count = 0
        protein_seq_count = 0
        chromosomal_location_count = 0
        locus_count = 0
        gene_seq_count = 0
        gene_description_count = 0

        for record in records:
            subject = record.get("subject", {})

            # protein function
            function = subject.get("function")
            if function:
                functions.append(function)

            # molecular weight
            mw = subject.get("molecular_weight")
            if mw is not None:
                molecular_weights.append(mw)

            if subject.get("transmembrane_region"):
                transmembrane_count += 1
            if subject.get("signal_region"):
                signal_region_count += 1
            if subject.get("protein_seq"):
                protein_seq_count += 1
            if subject.get("chromosomal_location"):
                chromosomal_location_count += 1
            if subject.get("locus"):
                locus_count += 1
            if subject.get("gene_seq"):
                gene_seq_count += 1
            if subject.get("gene_description"):
                gene_description_count += 1

        total_records = len(records)

        return {
            "function_stats": {
                "record_count": len(functions),
                "unique_count": len(set(functions)),
                "percentage": round((len(functions) / total_records) * 100, 2)
                if total_records > 0
                else 0,
            },
            "molecular_weight_stats": self._safe_get_numeric_stats(molecular_weights),
            "transmembrane_region_count": transmembrane_count,
            "signal_region_count": signal_region_count,
            "protein_seq_count": protein_seq_count,
            "chromosomal_location_count": chromosomal_location_count,
            "locus_count": locus_count,
            "gene_seq_count": gene_seq_count,
            "gene_description_stats": {
                "record_count": gene_description_count,
                "percentage": round((gene_description_count / total_records) * 100, 2)
                if total_records > 0
                else 0,
            },
            "cellular_component_stats": self._count_unique_items_in_list_field(
                records, ["subject", "cellular_component"]
            ),
        }

    def _analyze_protein_biological_process_specific(
        self, records: List[Dict[str, Any]]
    ) -> Dict[str, Any]:
        """Analyze protein-biological process specific statistics."""
        return {}  # No specific analysis needed beyond publication stats

    def save_stats_report(self, stats: Dict[str, Any]) -> str:
        """Save the statistics report to JSON file."""
        report_path = os.path.join(
            self.report_dir, f"record_stats_{datetime.now().strftime('%Y%m%d_%H%M%S')}.json"
        )
        with open(report_path, "w") as f:
            json.dump(stats, f, indent=2, sort_keys=True)

        print(f"Statistics report saved to: {report_path}")
        return report_path

    def run_full_analysis(self) -> Dict[str, Any]:
        """Run complete statistical analysis and save report."""
        print("Starting comprehensive HMDB record analysis...")

        stats = self.generate_record_stats()
        if stats:
            self.save_stats_report(stats)
            print("Analysis complete!")
        else:
            print("Analysis failed - no data found.")

        return stats

    def run_full_analysis_with_duplicates(self) -> Dict[str, Any]:
        """Run complete statistical analysis and export duplicate records."""
        print("Starting comprehensive HMDB record analysis with duplicate export...")

        stats = self.generate_record_stats()
        if stats:
            self.save_stats_report(stats)

        duplicate_data = self.export_duplicate_records()
        if stats and duplicate_data:
            print("Analysis and duplicate export complete!")
        elif stats:
            print("Analysis complete, but duplicate export failed.")
        else:
            print("Analysis failed - no data found.")

        return {"stats": stats, "duplicates": duplicate_data}


if __name__ == "__main__":
    reporter = HMDBRecordStatsReporter()
    full_report = reporter.run_full_analysis_with_duplicates()
