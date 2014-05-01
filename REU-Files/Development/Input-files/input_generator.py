#!/usr/bin/python

import optparse
import sys

parser = optparse.OptionParser()

# first pulse
parser.add_option('--rr1',
                  dest='rr1',# default='0.0d0',
                  help='pulse 1: delay in periods from the start of time')
parser.add_option('--n1up',
                  dest='n1up',# default='20',
                  help='pulse 1: number of periods for ramp-up ')
parser.add_option('--n1plat',
                  dest='n1plat',# default='0',
                  help='pulse 1: number of periods on plateau')
parser.add_option('--n1down',
                  dest='n1down',# default='20',
                  help='pulse 1: number of periods for ramp-down')
parser.add_option('--ww1',
                  dest='ww1',# default='0.5d0',
                  help='pulse 1: frequency in atomic units')
parser.add_option('--ee1',
                  dest='ee1',# default='0.05338d0',
                  help='pulse 1: amplitude in atomic units')
parser.add_option('--alph1',
                  dest='alph1',# default='1.0d0',
                  help='pulse 1: factor by which e1 is multiplied')
parser.add_option('--s1up',
                  dest='shape1up',# default='s',
                  help='pulse 1: character variable for shape of ramp-up')
parser.add_option('--s1down',
                  dest='shape1down',# default='s',
                  help='pulse 1: character variable for shape of ramp-dpwm')
parser.add_option('--cep1',
                  dest='cep1',# default='0.0d0',
                  help='pulse 1: carrier envelope phase relative to sin(ww1*t)')

(opts, args) = parser.parse_args()

# output the input file to stdout
# should be redirected to pulse.inp

inputs = opts.__dict__.items()

# print the reminding
first = True
sys.stdout.write(' &pulse1  ')
for i in range(len(inputs)):
    if inputs[i][1] != None:
        if(not first):
            sys.stdout.write(', ')
        else:
            first = False
        sys.stdout.write(inputs[i][0] +' = ' + inputs[i][1])
        
sys.stdout.write(' /\n')

sys.stdout.write(' &pulse2  ' + 'alph2=0.0d0' + ' /\n')
