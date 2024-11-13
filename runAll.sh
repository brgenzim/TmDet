#!/bin/bash

ROOT="/bigdisk/users/tusi/pdbtm"

for code in `cat ../data/ids`
do
	while [ "$(squeue | grep -c run.sh)" -gt 500 ]
	do
		echo "waiting..."
		sleep 30
	done
	LOG="$ROOT/log/"`echo $code | cut -c 2-3`
	sbatch -c 2 --mem 2G -e "$LOG/$code.e.log" -o "$LOG/$code.o.log" --chdir "$ROOT/TmDet" run.sh $code
done
