#!/usr/bin

import sys

for i in range(len(sys.argv)):
    if(sys.argv[i][0] == '-' and sys.argv[i][1] == 'w'):
        print "the wave is", sys.argv[i+1]

