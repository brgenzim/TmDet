#!/bin/bash

function whereisit() {
    file=`grep -r " $1(" ./TmdetLib/ ./PdbLib/ ./TmdetUtils/ | grep -v svn | grep -v '\.h:' | cut -d':' -f 1`
}

function getfunctions() {
    f=`grep -oPz "(?s) $1.*?\n}" $2 | grep -a -o '[= (].*[a-z](' | grep -v free | grep -v calloc | tr -d ' ' | tr -d '=' | tr -d '(' | tail -n +2`
}

function search() {
    whereisit $1
    getfunctions $1 $file
    for i in $f
    do
        echo $i
        #search $i
    done
}

search $1