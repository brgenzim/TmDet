#!/bin/bash

cd ../log
for i in ?? ; do cat $i/*.o.log; done > tmdet.log
grep ' no ' tmdet.log | cut -d' ' -f 1 > nottm.res
grep ' yes ' tmdet.log | cut -d' ' -f 1 > tm.res
grep -f nottm.res ../data/pdbtm3.0.res | grep ' yes' | cut -d' ' -f 1 > fn.res
grep -f tm.res ../data/pdbtm3.0.res | grep ' no' | cut -d' ' -f 1 > fp.res
echo -n "number of false negatives: "; wc -l fn.res
echo -n "number of false positives: "; wc -l fp.res

for i in `cat ../data/ids`
do 
    echo `grep $i ../data/pdbtm3.0.res`" - "`grep $i tmdet.log`
done > comp.res

cd ../TmDet
for i in `cat ../log/fp.res`
do
    ./bin/tmdet -r -s -c $i
done