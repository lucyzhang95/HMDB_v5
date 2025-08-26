"""Cache utilities for storing and retrieving processed data."""
import json
import os
import pickle

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, filename: str) -> None:
    """Save object to pickle file in cache directory."""
    filepath = os.path.join(CACHE_DIR, filename)
    with open(filepath, "wb") as f:
        pickle.dump(obj, f)


def load_pickle(filename: str):
    """Load object from pickle file in cache directory."""
    filepath = os.path.join(CACHE_DIR, filename)
    return pickle.load(open(filepath, "rb")) if os.path.exists(filepath) else None


def save_json(obj, filename: str) -> None:
    """Save object to JSON file in cache directory."""
    filepath = os.path.join(CACHE_DIR, filename)
    with open(filepath, "w") as f:
        json.dump(obj, f, indent=4)


def cache_exists(filename: str) -> bool:
    """Check if cache file exists."""
    filepath = os.path.join(CACHE_DIR, filename)
    return os.path.exists(filepath)


def get_cache_path(filename: str) -> str:
    """Get full path to cache file."""
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
    """List all files in cache directory."""
    if not os.path.exists(CACHE_DIR):
        return []
    return [f for f in os.listdir(CACHE_DIR) if os.path.isfile(os.path.join(CACHE_DIR, f))]
