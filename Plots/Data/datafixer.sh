#!/bin/bash

for dir in S-s__2-36-2__0.5338d-1__0375*; do
    cd $dir
    cat overlap.out | python ../data_parser.py > intermediate_fixed.out
    cd ..
done