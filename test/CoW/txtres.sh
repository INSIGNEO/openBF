#!/bin/bash

for g in *results
do
	cd $g
for f in *.last
do
	cp -- "$f" "${f%.last}.txt"
done

cd ..
done 

cd ..
