#!/bin/bash

# This is a shell script that automates the running of multiple jobs. 
# To get the script running on your account, you should only have to 
# change the environment variables (the ones in all capital letters). 
# If you want to change the behavior of the code you WILL have to edit 
# some of code in the "begin script section"

#####
# TEST NAME
#####

TEST_DESCRIPTION=1p11-1.688d-2


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
INP_FILE_GEN_FP=$CODE_DIR_FP/Input-files/input_generator_2pulse.py

# the full path to the tdse*.inp file you are using
NUMERICS_INPUT_FP=$CODE_DIR_FP/Input-files/tdse-w2.inp

#####
# TEST PARAMETERS: changes these to create tests with different inputs
#####

# these are the test parameters that allow you to loop through multiple tests. 
# You may want to change these depending on your goals
PARAMETER_1=(0.3200 0.2700 0.2800 0.2900 0.3000 0.3100 0.3300 0.3400 0.3500 0.3550 0.3600 0.3625 0.3650 0.3675 0.3700 0.3725 0.3750 0.3775 0.3800 0.3825 0.3850 0.3875 0.3900 0.3950 0.4000 0.4100 0.4200 0.4300 0.4300 0.4400 0.4500 0.4600)
PARAMETER_2=(045.0d0 090.0d0 135.0d0)

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
		# compute the number of cycles for the plateau from the ramp up/down
		let "plat = 0"

		mkdir $WORK_DIR_FP/$TEST_DESCRIPTION-$p1-ss20-cep$p2-del000-0p225
		cd $WORK_DIR_FP/$TEST_DESCRIPTION-$p1-ss20-cep$p2-del000-0p225
			
		# copy and create input files
		cp $source_dir/$numerics_filename ./tdse.inp 
		python2.6 $INP_FILE_GEN_FP --ee1=1.688d-2 --ww1=$p1 --x1up=20.0d0 --x1plat=0.0d0 --x1down=20.0d0 --s1up=\'s\' --s1down=\'s\' --cep1=0.0d0 --alph2=0.225d0 --cep2=$p2 > pulse.inp

		# run the code 
		$source_dir/$code_filename > $code_filename-$p1-ss20-cep$p2-del000-0p225.log

		# go up a dir and run the next test
		cd ..
	done
done 
