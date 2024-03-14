#!/bin/bash

if [ -f build ]
then
	rm -rf build
fi

mkdir build
cd build
cmake ..
make && make install
cd ..

