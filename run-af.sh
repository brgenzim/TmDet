#!/bin/bash

AF_ROOT="/bigdisk/databases/AlphaFold/data"
TMDET_ROOT="/bigdisk/users/tusi/pdbtm/data/AlphaFold"

echo "./bin/tmdet -r -pi $AF_ROOT/$1/AF-$2-F1-model_v2.cif.gz -po $TMDET_ROOT/$1/$2.cif.gz -xo $TMDET_ROOT/$1/$2.xml -fr -minht"
./bin/tmdet -r -pi $AF_ROOT/$1/AF-$2-F1-model_v2.cif.gz -po $TMDET_ROOT/$1/$2.cif.gz -xo $TMDET_ROOT/$1/$2.xml -fr -minht 10
