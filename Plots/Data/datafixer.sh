#!/bin/bash

for dir in S-s*; do
    cd $dir
    cat overlap.out | python ../data_parser.py > intermediate_fixed.out
    cd ..
done