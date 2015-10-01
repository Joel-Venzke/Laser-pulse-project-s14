#!/bin/bash

# This is a shell script that automates the running of multiple jobs. 
# To get the script running on your account, you should only have to 
# change the environment variables (the ones in all capital letters). 
# If you want to change the behavior of the code you WILL have to edit 
# some of code in the "begin script section"

#####
# TEST NAME
#####

TEST_DESCRIPTION=3pulse


#####
# SET THE ENVIRONMENT VARIABLES: you should only have to set these once. 
# The information below is for John Emmons' account, so make sure you 
# change the code to match where you have put the code
#####

# the full path to your Development directory
# to find the path, type "pwd" in your Development dir 
CODE_DIR_FP=/home1/02971/jvenzke/Laser-pulse-project-s14/REU-Files/Development

# the full path to your Work direcoty
# to find the path, type "cdw" to go to your work directory, then "pwd" 
WORK_DIR_FP=/work/02971/jvenzke

# the full path to the compiled code 
COMPILED_CODE_FP=$CODE_DIR_FP/work-29

# the full path to the input file generator
INP_FILE_GEN_FP=$CODE_DIR_FP/Input-files/input_generator_3pulse.py

# the full path to the tdse*.inp file you are using
NUMERICS_INPUT_FP=$CODE_DIR_FP/Input-files/tdse-3pulse-IR-klaus.inp
NE_INPUT_FP=$CODE_DIR_FP/Input-files/ne.wfn

#####
# TEST PARAMETERS: changes these to create tests with different inputs
#####

# these are the test parameters that allow you to loop through multiple tests. 
# You may want to change these depending on your goals
#PARAMETER_1=( 0.005d0 0.010d0 0.015d0 0.020d0 )
#PARAMETER_1=( 0.0215d0 0.0212d0 0.0217d0 )
PARAMETER_1=( 0.02155d0 )

#####
# BEGIN THE SCRIPT BELOW: you will have to edit the code below if you want to 
# change the behavior of the code or change the way tests are run.
#####

# get the name of the compiled code and numerics file
code_filename=$(basename $COMPILED_CODE_FP)
numerics_filename=$(basename $NUMERICS_INPUT_FP) 
ne_filename=$(basename $NE_INPUT_FP) 

# copy code to the work dir
mkdir $WORK_DIR_FP/$TEST_DESCRIPTION-src

# get the full path to the source (src) directory that just was created
# you should NOT change this code!
cd $WORK_DIR_FP/$TEST_DESCRIPTION-src   
source_dir=$(pwd)

cp $COMPILED_CODE_FP .
cp $NUMERICS_INPUT_FP . 
cp $NE_INPUT_FP .

# go up a directory and run the tests
cd ..

# begin the main loop for the tests. This will do all possible combinations of the elements of the parameter lists
for p1 in ${PARAMETER_1[*]}; do
	# compute the number of cycles for the plateau from the ramp up/down
	let "plat = 0"
	TMP_FILENAME=$TEST_DESCRIPTION-0.0d0-$p1
	mkdir $WORK_DIR_FP/${TMP_FILENAME}
	cd $WORK_DIR_FP/${TMP_FILENAME}
		
	# copy and create input files
	cp $source_dir/$numerics_filename ./tdse.inp 
	cp $source_dir/$ne_filename ./ne.wfn 

	# create pulse.inp
	python2.6 $INP_FILE_GEN_FP --alph1=0.0d0 --ee1=5.338d-3 --ww1=0.7092d0 --x1up=125.0d0 --x1plat=0.0d0 --x1down=125.0d0 --s1up=\'s\' --s1down=\'s\' --cep1=0.0d0 \
	--alph2=$p1 --rr2=0.0d0 \
	--alph3=0.d0 --ee3=5.338d-4 --ww3=0.06d0 --x3up=3.00d0 --x3plat=0.0d0 --x3down=3.00d0 > pulse.inp

	# run the code 
	$source_dir/$code_filename > $code_filename-${TMP_FILENAME}.log

	# go up a dir and run the next test
	cd ..
done 
