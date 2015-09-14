#!/bin/bash
#SBATCH -J 3pulse_rr2.nam
#SBATCH -o 3pulse_rr2.log
#SBATCH -e 3pulse_rr2.err 
#SBATCH -N  1  -n  1
#SBATCH -p  normal
#SBATCH -t  48:00:00

export OMP_NUM_THREADS=16
export KMP_STACKSIZE=16m

# you need to cd into the Job-scripts directory or wherever you have put the corresponding job shell script
cd /home1/02971/jvenzke/Laser-pulse-project-s14/REU-Files/Job-scripts

# This is a shell script that automates the running of multiple jobs. 
# To get the script running on your account, you should only have to 
# change the environment variables (the ones in all capital letters). 
# If you want to change the behavior of the code you WILL have to edit 
# some of code in the "begin script section"

#####
# TEST NAME
#####

TEST_DESCRIPTION=3p


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
WORK_DIR_FP=/work/02971/jvenzke/3pulse

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
PARAMETER_1=( 0.685d0 )
PARAMETER_2=( 0.00d0  0.05d0  0.10d0  0.15d0  0.20d0  0.25d0  0.30d0  0.35d0  0.40d0  0.45d0  0.50d0  0.55d0  0.60d0  0.65d0  0.70d0  0.75d0  0.80d0  0.85d0  0.90d0  0.95d0)

#####
# BEGIN THE SCRIPT BELOW: you will have to edit the code below if you want to 
# change the behavior of the code or change the way tests are run.
#####

# get the name of the compiled code and numerics file
code_filename=$(basename $COMPILED_CODE_FP)
numerics_filename=$(basename $NUMERICS_INPUT_FP) 
ne_filename=$(basename $NE_INPUT_FP) 

# copy code to the work dir
mkdir -p $WORK_DIR_FP/$TEST_DESCRIPTION-src

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
	for p2 in ${PARAMETER_2[*]}; do
		# compute the number of cycles for the plateau from the ramp up/down
		let "plat = 0"
		echo TMP_FILENAME=$TEST_DESCRIPTION-$p1-$p2
		mkdir -p $WORK_DIR_FP/${TMP_FILENAME}
		cd $WORK_DIR_FP/${TMP_FILENAME}
			
		# copy and create input files
		cp $source_dir/$numerics_filename ./tdse.inp 
		cp $source_dir/$ne_filename ./ne.wfn 

		# create pulse.inp
		python2.6 $INP_FILE_GEN_FP --alph1=1.0d0 --ee1=5.338d-3 --ww1=$p1 --x1up=125.0d0 --x1plat=0.0d0 --x1down=125.0d0 --s1up=\'s\' --s1down=\'s\' --cep1=0.0d0 \
		--alph2=0.02155d0 --rr2=0.0d0 \
		--alph3=0.d0 --ee3=5.338d-4 --ww3=0.06d0 --x3up=3.00d0 --x3plat=0.0d0 --x3down=3.00d0 > pulse.inp

		# run the code 
		$source_dir/$code_filename > run.log

		# go up a dir and run the next test
		cd ..
	done
done 
