set -e
if [[ "$#" == 0 ]]; then
  {
    echo "Usage: bash $0 <ccd-data-directory-path>"
    echo "Desitination directory must be passed as argument to download CCD data"
  }  > /dev/stderr
  exit 1
fi
DEST_DIR="$1"

mkdir -p "$DEST_DIR"
curl --fail -Lo "$DEST_DIR/components.cif.gz" \
  "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"

