# dump.sh: download all URLs passed as arguments
set -euo pipefail

urls=(
  "https://hmdb.ca/system/downloads/current/hmdb_metabolites.zip"
  "https://hmdb.ca/system/downloads/current/hmdb_proteins.zip"
  "https://smpdb.ca/downloads/smpdb_pathways.csv.zip"
)

for url in "${urls[@]}"; do
  echo "Downloading $url …"
  wget --continue --timestamping --progress=bar:force "$url"

  echo
done

echo "All HMDB downloads are complete."

# Usage:
# change the directory to where the script is located
# chmod +x dump.sh
# ./dump.sh