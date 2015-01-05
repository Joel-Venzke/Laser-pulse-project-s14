#!/bin/bash

# This is a shell script that automates the running of multiple jobs. 
# To get the script running on your account, you should only have to 
# change the environment variables (the ones in all capital letters). 
# If you want to change the behavior of the code you WILL have to edit 
# some of code in the "begin script section"

#####
# TEST NAME
#####

TEST_DESCRIPTION=1.750d-1-joel_test_0350_short_pulse


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
COMPILED_CODE_FP=$CODE_DIR_FP/work25

# the full path to the input file generator
INP_FILE_GEN_FP=$CODE_DIR_FP/Input-files/input_generator.py

# the full path to the tdse*.inp file you are using
NUMERICS_INPUT_FP=$CODE_DIR_FP/Input-files/tdse-3.5cycle.inp

#####
# TEST PARAMETERS: changes these to create tests with different inputs
#####

# these are the test parameters that allow you to loop through multiple tests. 
# You may want to change these depending on your goals
PARAMETER_1=(1.0d0 1.5d0 2.0d0 2.5d0 3.0d0 3.5d0)
PARAMETER_2=(g)
PARAMETER_3=(90.0d0)

#####
# BEGIN THE SCRIPT BELOW: you will have to edit the code below if you want to 
# change the behavior of the code or change the way tests are run.
#####

# get the name of the compiled code and numerics file
code_filename=$(basename $COMPILED_CODE_FP)
numerics_filename=$(basename $NUMERICS_INPUT_FP) 

# copy code to the work dir
mkdir $WORK_DIR_FP/$TEST_DESCRIPTION-src

# get the full path to the source (src) directory that just was created
# you should NOT change this code!
cd $WORK_DIR_FP/$TEST_DESCRIPTION-src   
source_dir=$(pwd)

cp $COMPILED_CODE_FP .
cp $NUMERICS_INPUT_FP . 

# go up a directory and run the tests
cd ..

# begin the main loop for the tests. This will do all possible combinations of the elements of the parameter lists
for p1 in ${PARAMETER_1[*]}; do
    for p2 in ${PARAMETER_2[*]}; do	
    	for p3 in ${PARAMETERS_3[*]}; do

# compute the number of cycles for the plateau from the ramp up/down
		let "plat = 0"

		mkdir $WORK_DIR_FP/$TEST_DESCRIPTION_CEP-$p3-$p1-$plat-$p1-$p2-$p2
		cd $WORK_DIR_FP/$TEST_DESCRIPTION_CEP-$p3-$p1-$plat-$p1-$p2-$p2
		
# copy and create input files
		cp $source_dir/$numerics_filename ./tdse.inp 
		python2.6 $INP_FILE_GEN_FP --ee1=1.750d-1 --ww1=0.350d0 --x1up=$p1 --x1plat=$plat --x1down=$p1 --s1up=\'$p2\' --s1down=\'$p2\' --cep1=$p3 > pulse.inp

# run the code 
		$source_dir/$code_filename > $code_filename_CEP-$p3-$p1-$plat-$p1-$p2-$p2.log 

# go up a dir and run the next test
		cd ..
		done
    done
done 
