# Start debug server: Mol* visualization
# usage: bash tests/bin/dbgsrv.sh PDB_CODE RADIUS THICKNESS CIF_FILE

php -S localhost:8081 -t tests/web \
    -d pdbCode="$1" \
    -d radius="$2" \
    -d thickness="$3" \
    -d cifFile="$(readlink -f $4)"
