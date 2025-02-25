#!/bin/bash

AF_ROOT="/bigdisk/users/tusi/pdbtm/data/AlphaFold/new/test/afdb"
TMDET_ROOT="/bigdisk/users/tusi/pdbtm/data/AlphaFold/new/test/tmdet"

echo "./bin/tmdet -r -pi $AF_ROOT/AF-$1-F1-model_v4.cif -po $TMDET_ROOT/$1.cif.gz -xo $TMDET_ROOT/$1.xml -fr -minht 13"
./bin/tmdet -r -pi $AF_ROOT/AF-$1-F1-model_v4.cif -po $TMDET_ROOT/$1.cif.gz -xo $TMDET_ROOT/$1.xml -fr -minht 13
