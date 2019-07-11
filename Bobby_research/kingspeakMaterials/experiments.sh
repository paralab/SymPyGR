#!/bin/bash

for filename in ~/research/rhs_codes/*.cpp; do
	echo $filename
	cp $filename ../src/rhs.cpp
	make all -j4
	./rhsTest 0 4 50 0 > "./results/$(basename "$filename" .cpp).txt"
done
