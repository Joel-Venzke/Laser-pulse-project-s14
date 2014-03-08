#/bin/bash


TEST_DESCRIPTION=Env-test
COMPILED_CODE=work15b-all


PARAMETER_1=(s g t)
PARAMETER_2=(s g t)
NUMERICS_INPUT=tdse-40cycle.inp

CODE_DIR=/home1/02603/jemmons/physics-capstone/REU-Files/Development
WORK_DIR=/work/02603/jemmons

# copy code to the work dir
mkdir $WORK_DIR/$TEST_DESCRIPTION-src
cd $WORK_DIR/$TEST_DESCRIPTION-src
local_code_path=$(pwd)

cp $CODE_DIR/$COMPILED_CODE .

# go up a dir and run the tests
cd ..

for p1 in ${PARAMETER_1[*]}; do
    for p2 in ${PARAMETER_2[*]}; do	

	mkdir $WORK_DIR/$TEST_DESCRIPTION-$p1-$p2
	cd $WORK_DIR/$TEST_DESCRIPTION-$p1-$p2
		
# copy and create input files
	cp $CODE_DIR/$NUMERICS_INPUT ./tdse.inp 
	python $CODE_DIR/input_generator.py --ee1=1.688d-1 --ww1=0.4d0 --n1up=2 --n1down=2 --shape1up=$p1 --shape1down=$p2 > pulse.inp

# run the code 
	$local_code_path/$COMPILED_CODE > $COMPILED_CODE.log 

# go up a dir and run the next test
	cd ..

    done
done 