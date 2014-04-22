#!/usr/bin

import sys

# read input from stdin

file_data = sys.stdin.readlines()

#print file_data

j = 100249
for i in range(len(file_data)):
    if file_data[i][0] == '*':
        line = file_data[i].split('\t')
        line[0] = str(j)
        j += 250
        sys.stdout.write('\t'.join(line))
    else:
        sys.stdout.write(file_data[i])
