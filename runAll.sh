#!/bin/bash

ROOT="/bigdisk/users/tusi/pdbtm"

while IFS= read -r line
do
	while [ "$(squeue | grep -c run.sh)" -gt 500 ]
	do
		echo "waiting..."
		sleep 30
	done
	code=`echo $line | cut -c 1-4`
	args=`echo $line | cut -c 5-`
	LOG="$ROOT/log/"`echo $code | cut -c 2-3`
	sbatch -p rall -c 2 --mem 8G -e "$LOG/$code.e.log" -o "$LOG/$code.o.log" --chdir "$ROOT/TmDet" run.sh $code $args
done  < $1
