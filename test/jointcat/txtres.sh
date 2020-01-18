#!/bin/bash

cd conjunction_results

for f in *.last
do
	cp -- "$f" "${f%.last}.txt"
done

cd ..
