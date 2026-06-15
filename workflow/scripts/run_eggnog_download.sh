#!/bin/bash
set -euo pipefail
# Patch BASE_URL in the installed download script (eggnogdb.embl.de is dead).
DL=$(which download_eggnog_data.py)
sed -i 's|http://eggnogdb.embl.de|http://eggnog5.embl.de|g' "$DL"
# Now run the patched downloader
mkdir -p resources/eggnog_data
download_eggnog_data.py --data_dir resources/eggnog_data/ -y
