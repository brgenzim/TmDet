#!/bin/bash

ROOT="/bigdisk/users/tusi/pdbtm"
AF_ROOT="/bigdisk/databases/AlphaFold/data"
TMDET_ROOT="/bigdisk/users/tusi/pdbtm/data/AlphaFold"

mkdir -p $TMDET_ROOT/$1

for code in `cat $AF_ROOT/$1/ids.txt`
do
	while [ "$(squeue | grep -c run.sh)" -gt 500 ]
	do
		echo "waiting..."
		sleep 30
	done
	LOG="$ROOT/log/"`echo $code | cut -c 2-3`
	echo -n "$code "
	sbatch -p rall -c 2 --mem 8G -e "$LOG/$code.e.log" -o "$LOG/$code.o.log" --chdir "$ROOT/TmDet" run-af.sh $1 $code
done
