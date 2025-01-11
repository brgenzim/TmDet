#!/bin/bash

ROOT="/bigdisk/users/tusi/pdbtm"

while IFS= read -r line
do
	while [ "$(squeue | grep -c run.sh)" -gt 500 ]
	do
		echo "waiting..."
		sleep 30
	done
	code=`echo $line | cut -d' ' -f 1`
	chains=`echo $line | cut -d' ' -f 2`
	LOG="$ROOT/log/"`echo $code | cut -c 2-3`
	#sbatch -p btwins -c 2 --mem 2G -e "$LOG/$code.e.log" -o "$LOG/$code.o.log" --chdir "$ROOT/TmDet" 
	echo $code
	./run.sh $code $chains
done  < $1
