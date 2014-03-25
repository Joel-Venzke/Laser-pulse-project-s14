#!/bin/bash

TEST_DESCRIPTION=Env-test-0350
COMPILED_CODE=work15b-all

PARAMETER_1=(s t)
PARAMETER_2=(s t)
NUMERICS_INPUT=tdse-40cycle.inp

CODE_DIR=/home1/02869/smbuczek/physics-capstone/REU-Files/Development
WORK_DIR=/work/02869/smbuczek

# copy code to the work dir
mkdir $WORK_DIR/$TEST_DESCRIPTION-src
cd $WORK_DIR/$TEST_DESCRIPTION-src
local_code_path=$(pwd)

cp $CODE_DIR/$COMPILED_CODE .
cp $CODE_DIR/$NUMERICS_INPUT .

# go up a dir and run the tests
cd ..

for p1 in ${PARAMETER_1[*]}; do
    for p2 in ${PARAMETER_2[*]}; do	

	mkdir $WORK_DIR/$TEST_DESCRIPTION-$p1-$p2
	cd $WORK_DIR/$TEST_DESCRIPTION-$p1-$p2
		
# copy and create input files
	cp $local_code_path/$NUMERICS_INPUT ./tdse.inp 
	python $CODE_DIR/input_generator.py --ee1=1.0676d-1 --ww1=0.350d0 --n1up=2 --n1plat=36 --n1down=2 --s1up=\'$p1\' --s1down=\'$p2\' > pulse.inp

# run the code 
	$local_code_path/$COMPILED_CODE > $COMPILED_CODE-$p1-$p2.log 

# go up a dir and run the next test
	cd ..

    done
done 