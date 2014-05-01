TARGET = work15b-all
COMPILER = ifort

# production/debugging compiler flags
# only one should ever be uncommented at a time
CFLAGS = -openmp
#CFLAGS = -C -traceback 


all:
	$(COMPILER) -o $(TARGET) $(CFLAGS) $(TARGET).f
	rm -rf *.mod *.pyc

clean:
	rm -rf work15b-all *.mod read.me *~ *.out