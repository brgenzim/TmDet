#!/bin/bash

if [ "$1" != "$2" ]; then
    chains=" -uc $2"
else
    chains=""
fi

./bin/tmdet -r -c $1 $chains
