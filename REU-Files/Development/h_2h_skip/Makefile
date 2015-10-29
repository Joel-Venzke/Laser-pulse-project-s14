TARGET = work-29-ne8-rr2-fixed
COMPILER = gfortran

# production/debugging compiler flags
# only one should ever be uncommented at a time
# CFLAGS = -openmp
CFLAGS = 


all:
	$(COMPILER) -o $(TARGET) $(CFLAGS) $(TARGET).f
	rm -rf *.mod *.pyc

clean:
	rm -rf run/*.mod run/read.me run/*~ run/*.out