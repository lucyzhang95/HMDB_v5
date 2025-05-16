import requests
from bs4 import BeautifulSoup
import re


def get_release(self):
    """Retrieves the current version number of the HMDB from the downloads page.

    :returns: The current version number of the HMDB if found, else None.
    """
    url = "https://hmdb.ca/downloads"

    try:
        response = requests.get(url, timeout=10)
        response.raise_for_status()
    except requests.RequestException as e:
        print(f"Failed to fetch HMDB page: {e}")
        return None

    soup = BeautifulSoup(response.text, "html.parser")
    active_li = soup.select_one("#download_tabs li.active a")
    if active_li:
        match = re.search(r"\(([\d.]+)\)", active_li.text)
        if match:
            current_version = match.group(1)
            return current_version

    return None

if __name__ == "__main__":
    version = get_release(self)
    print("Current HMDB version: {}".format(version))

