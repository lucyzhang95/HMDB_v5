"""Cache utilities for storing and retrieving processed data."""
import json
import os
import pickle

CACHE_DIR = os.path.join(os.getcwd(), "cache")
REC_DIR = os.path.join(os.getcwd(), "records")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, filename: str) -> None:
    """Save object to pickle file in cache directory."""
    filepath = os.path.join(CACHE_DIR, filename)
    with open(filepath, "wb") as in_f:
        pickle.dump(obj, in_f)


def load_pickle(filename: str):
    """Load object from pickle file in cache directory."""
    filepath = os.path.join(CACHE_DIR, filename)
    return pickle.load(open(filepath, "rb")) if os.path.exists(filepath) else None


def save_json(obj, filename: str) -> None:
    """Save object to JSON file in cache directory."""
    filepath = os.path.join(REC_DIR, filename)
    with open(filepath, "w") as out_f:
        json.dump(obj, out_f, indent=4)


def save_jsonl(records: List[Dict[str, Any]], f_name: str):
    """Saves records to a JSONL file with standardized formatting."""
    if not records:
        print("!!! Warning: No records provided for JSONL export.")
        return None

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    base_name = f_name.replace(".jsonl", "") if f_name.endswith(".jsonl") else f_name
    jsonl_filename = f"{base_name}_{timestamp}.jsonl"
    jsonl_path = os.path.join(REC_DIR, jsonl_filename)

    exported_count = 0
    with open(jsonl_path, "w", encoding="utf-8") as f:
        for record in records:
            clean_record = _standardize_record(record)
            json.dump(clean_record, f, ensure_ascii=False, separators=(",", ":"))
            f.write("\n")
            exported_count += 1

    print(f"-> Exported {exported_count} records to {jsonl_path}")
    return jsonl_path


def load_jsonl(f_name: str) -> List[Dict[str, Any]]:
    """Load records from JSONL file."""
    path = os.path.join(REC_DIR, f_name)

    if not os.path.exists(path):
        return []

    records = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    records.append(json.loads(line))
                except json.JSONDecodeError:
                    continue

    print(f"-> Loaded {len(records)} records from {f_name}")
    return records


def _standardize_record(record: Dict[str, Any]) -> Dict[str, Any]:
    """Standardizes a record by cleaning up empty values and formatting xrefs into a list."""
    clean_record = record.copy()

    for entity_key in ["subject", "object"]:
        if entity_key in clean_record and isinstance(clean_record[entity_key], dict):
            entity = clean_record[entity_key].copy()
            if "xrefs" in entity and isinstance(entity["xrefs"], dict):
                entity["xrefs"] = [v for k, v in entity["xrefs"].items() if v]
            clean_record[entity_key] = entity

    return clean_record


def cache_exists(filename: str) -> bool:
    """Check if the cache file exists."""
    filepath = os.path.join(CACHE_DIR, filename)
    return os.path.exists(filepath)


def get_cache_path(filename: str) -> str:
    """Get the full path to cache the file."""
    return os.path.join(CACHE_DIR, filename)


def clear_cache() -> None:
    """Remove all cache files."""
    import shutil

    if os.path.exists(CACHE_DIR):
        shutil.rmtree(CACHE_DIR)
        os.makedirs(CACHE_DIR, exist_ok=True)


def get_cache_size() -> dict:
    """Get size information about cache files."""
    cache_info = {}
    if not os.path.exists(CACHE_DIR):
        return cache_info

    for filename in os.listdir(CACHE_DIR):
        filepath = os.path.join(CACHE_DIR, filename)
        if os.path.isfile(filepath):
            size = os.path.getsize(filepath)
            cache_info[filename] = {"size_bytes": size, "size_mb": round(size / (1024 * 1024), 2)}

    return cache_info


def list_cache_files() -> list:
    """List all files in the cache directory."""
    if not os.path.exists(CACHE_DIR):
        return []
    return [f for f in os.listdir(CACHE_DIR) if os.path.isfile(os.path.join(CACHE_DIR, f))]
