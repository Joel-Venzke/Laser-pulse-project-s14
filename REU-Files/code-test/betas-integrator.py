#!/usr/bin
import sys

FILE_NAME = "betas.out"

inFile = open(FILE_NAME,'r')
file_string = inFile.read()

data_all = file_string.split('\n')

data = []
for line in data_all:
    pair = line.split('  ')[1:3]

    # append the element if it is nonempty
    if(bool(pair)):
        pair = [float(pair[0]), float( pair[1])]
        data.append(pair)

# linearly extrapolate the first point (ie. 0)
x1 = data[0][0]
y1 = data[0][1]
x2 = data[1][0]
y2 = data[1][1]

y0 = y1 - ((y2 - y1)/(x2 - x1)) * x1

# insert 0 at the front of the list
data.insert(0, [0.0, y0])


# debugging output
#for pair in data:
#    sys.stdout.write(str(pair[0]) + " " + str(pair[1]) + "\n")


# compute the integral using trapezoids
integral = 0.0;
for i in range(len(data) - 1):
    integral += ((data[i][1] + data[i + 1][1]) / 2) * (data[i + 1][0] - data[i][0])

print "The integral is", integral
