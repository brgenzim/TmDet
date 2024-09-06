set -e
if [[ "$#" == 0 ]]; then
  {
    echo "Usage: bash $0 <ccd-data-directory-path>"
    echo "Desitination directory must be passed as argument to download CCD data"
  }  > /dev/stderr
  exit 1
fi
DEST_DIR="$1"
COMPONENTS="$DEST_DIR/components.cif.gz"

mkdir -p "$DEST_DIR"
curl --fail -Lo "$COMPONENTS" \
  "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"

cat <<EOF

Use this command line to split components file:
   fragment_cif -i $COMPONENTS -d $DEST_DIR -s &> fragment.log

EOF
