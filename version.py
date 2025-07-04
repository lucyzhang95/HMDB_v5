import re
import time

import requests
from bs4 import BeautifulSoup


def get_release(max_retries=3, base_delay=2):
    """Retrieves the current version number of the HMDB from the downloads page.
    HMDB website is very slow and unstable sometimes,
    so need to wait for at least 30 sec.

    :returns: The current version number of the HMDB if found, else None.
    """
    url = "https://hmdb.ca/downloads"
    headers = {
        "User-Agent": (
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
            "AppleWebKit/537.36 (KHTML, like Gecko) "
            "Chrome/124.0 Safari/537.36"
        )
    }
    for attempt in range(1, max_retries + 1):
        try:
            print(f"Attempt {attempt} to fetch HMDB version...")
            response = requests.get(url, headers=headers, timeout=30)
            response.raise_for_status()

            soup = BeautifulSoup(response.text, "html.parser")
            active_li = soup.select_one("#download_tabs li.active a")
            if active_li:
                match = re.search(r"\(([\d.]+)\)", active_li.text)
                if match:
                    version = match.group(1)
                    return version
            print("HMDB version pattern not found.")
            return None

        except requests.RequestException as e:
            print(
                f"HMDB attempt {attempt} failed ({e}); retrying in {base_delay * (2 ** (attempt - 1))}s"
            )
            time.sleep(base_delay * (2 ** (attempt - 1)))

    print("Failed to fetch HMDB version after multiple attempts.")
    return None


if __name__ == "__main__":
    version = get_release()
    print("Current HMDB version", version)
