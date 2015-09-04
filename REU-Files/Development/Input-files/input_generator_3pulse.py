#!/usr/bin/python

import optparse
import sys

parser = optparse.OptionParser()

# first pulse
parser.add_option('--rr1',
                  dest='rr1',# default='0.0d0',
                  help='pulse 1: delay in periods from the start of time')
parser.add_option('--x1up',
                  dest='x1up',# default='20',
                  help='pulse 1: number of periods for ramp-up ')
parser.add_option('--x1plat',
                  dest='x1plat',# default='0',
                  help='pulse 1: number of periods on plateau')
parser.add_option('--x1down',
                  dest='x1down',# default='20',
                  help='pulse 1: number of periods for ramp-down')
parser.add_option('--ww1',
                  dest='ww1',# default='0.5d0',
                  help='pulse 1: frequency in atomic units')
parser.add_option('--ee1',
                  dest='ee1',# default='0.05338d0',
                  help='pulse 1: amplitude in atomic units')
parser.add_option('--alph1',
                  dest='alph1',# default='1.0d0',
                  help='pulse 1: factor by which ee1 is multiplied')
parser.add_option('--s1up',
                  dest='shape1up',# default='s',
                  help='pulse 1: character variable for shape of ramp-up')
parser.add_option('--s1down',
                  dest='shape1down',# default='s',
                  help='pulse 1: character variable for shape of ramp-dpwm')
parser.add_option('--cep1',
                  dest='cep1',# default='0.0d0',
                  help='pulse 1: carrier envelope phase relative to sin(ww1*t)')

# second pulse
parser.add_option('--rr2',
                  dest='rr2',# default='0.0d0',
                  help='pulse 2: delay in periods from the start of time')
parser.add_option('--x2up',
                  dest='x2up',# default='20',
                  help='pulse 2: number of periods for ramp-up ')
parser.add_option('--x2plat',
                  dest='x2plat',# default='0',
                  help='pulse 2: number of periods on plateau')
parser.add_option('--x2down',
                  dest='x2down',# default='20',
                  help='pulse 2: number of periods for ramp-down')
parser.add_option('--ww2',
                  dest='ww2',# default='0.5d0',
                  help='pulse 2: frequency in atomic units')
parser.add_option('--ee2',
                  dest='ee2',# default='0.05338d0',
                  help='pulse 2: amplitude in atomic units')
parser.add_option('--alph2',
                  dest='alph2', default='0.0d0',
                  help='pulse 2: factor by which ee2 is multiplied')
parser.add_option('--s2up',
                  dest='shape2up',# default='s',
                  help='pulse 2: character variable for shape of ramp-up')
parser.add_option('--s2down',
                  dest='shape2down',# default='s',
                  help='pulse 2: character variable for shape of ramp-dpwm')
parser.add_option('--cep2',
                  dest='cep2',# default='0.0d0',
                  help='pulse 2: carrier envelope phase relative to sin(ww1*t)')


# third pulse
parser.add_option('--rr3',
                  dest='rr3',# default='0.0d0',
                  help='pulse 3: delay in periods from the start of time')
parser.add_option('--x3up',
                  dest='x3up',# default='20',
                  help='pulse 3: number of periods for ramp-up ')
parser.add_option('--x3plat',
                  dest='x3plat',# default='0',
                  help='pulse 3: number of periods on plateau')
parser.add_option('--x3down',
                  dest='x3down',# default='20',
                  help='pulse 3: number of periods for ramp-down')
parser.add_option('--ww3',
                  dest='ww3',# default='0.5d0',
                  help='pulse 3: frequency in atomic units')
parser.add_option('--ee3',
                  dest='ee3',# default='0.05338d0',
                  help='pulse 3: amplitude in atomic units')
parser.add_option('--alph3',
                  dest='alph3', default='0.0d0',
                  help='pulse 3: factor by which ee3 is multiplied')
parser.add_option('--s3up',
                  dest='shape3up',# default='s',
                  help='pulse 3: character variable for shape of ramp-up')
parser.add_option('--s3down',
                  dest='shape3down',# default='s',
                  help='pulse 3: character variable for shape of ramp-dpwm')
parser.add_option('--cep3',
                  dest='cep3',# default='0.0d0',
                  help='pulse 3: carrier envelope phase relative to sin(ww1*t)')

(opts, args) = parser.parse_args()

# output the input file to stdout
# should be redirected to pulse.inp

inputs = opts.__dict__.items()

def check_pulse_1(input_val):
    pulse1 = ["rr1", "x1up", "x1plat", "x1down", "ww1", "ee1", "alph1", "shape1up", "shape1down", "cep1"]
    if input_val in pulse1:
        return True
    return False

def check_pulse_2(input_val):
    pulse2 = ["rr2", "x2up", "x2plat", "x2down", "ww2", "ee2", "alph2", "shape2up", "shape2down", "cep2"]
    if input_val in pulse2:
        return True
    return False

def check_pulse_3(input_val):
    pulse3 = ["rr3", "x3up", "x3plat", "x3down", "ww3", "ee3", "alph3", "shape3up", "shape3down", "cep3"]
    if input_val in pulse3:
        return True
    return False


# print pulse 1
first = True
sys.stdout.write(' &pulse1  ')
for i in range(len(inputs)):
    if inputs[i][1] != None:
        if (check_pulse_1(inputs[i][0])):
            if(not first):
                sys.stdout.write(', ')
            else:
                first = False
            sys.stdout.write(inputs[i][0] +' = ' + inputs[i][1])
        
sys.stdout.write(' /\n')

# print pulse 2
first = True
sys.stdout.write(' &pulse2  ')
for i in range(len(inputs)):
    if inputs[i][1] != None:
        if (check_pulse_2(inputs[i][0])):
            if(not first):
                sys.stdout.write(', ')
            else:
                first = False
            sys.stdout.write(inputs[i][0] +' = ' + inputs[i][1])
sys.stdout.write(' /\n')

# print pulse 3
first = True
sys.stdout.write(' &pulse3  ')
for i in range(len(inputs)):
    if inputs[i][1] != None:
        if (check_pulse_3(inputs[i][0])):
            if(not first):
                sys.stdout.write(', ')
            else:
                first = False
            sys.stdout.write(inputs[i][0] +' = ' + inputs[i][1])
sys.stdout.write(' /\n')
