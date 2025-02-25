#!/bin/bash

ROOT="/bigdisk/users/tusi/pdbtm"

for code in `cat $1`
do
	while [ "$(squeue | grep -c run-af2.sh)" -gt 500 ]
	do
		echo "waiting..."
		sleep 30
	done
	LOG="$ROOT/log/af/"`echo $code | cut -c 2-5`
	echo -n "$code "
	sbatch -p rall -c 2 --mem 8G -e "$LOG/$code.e.log" -o "$LOG/$code.o.log" --chdir "$ROOT/TmDet" run-af2.sh $code
done
