READ ME

The files you will need to run the jobs are in this folder.
The files named titan and stampeded are set up to run on titan and stampeded respectively. 
You will need the input_generator_2pulse.py file for both computers.
In the uni files, you will need to update the wall time and the name of the files you wish to have for .err .log and the executable (ending in .sh)
The .uni file will also need the cd to be updates to where the .sh file lives
The .sh files will need changes to run on your setup. (explained bellow)
Both the titan and stampede examples are set up to run the latest job I sent data for.

Stampede.sh:
	You will need to update the following variables to here your code and other directories are is

	CODE_DIR_FP 
		where you keep your complied code
	WORK_DIR_FP
		Where you want your output folders to be made
	COMPILED_CODE_FP
		The version of code you want to run
	INP_FILE_GEN_FP
		Path to the input_generator_2pulse.py file
	NUMERICS_INPUT_FP
		path to the tdse.inp file you would like to use for this run
	PARAMETER_1=(0.360 0.375 0.390)
	PARAMETER_2=(000.0 045.0 090.0 135.0)
	PARAMETER_3=(0.0 0.5 1.0 1.5 2.0)
		the PARAMETER variables hold the list of variables you wish to test.
	now you can skip to the for loops
		each for loop will go the the respective PARAMETER array and store it in $p1, $p2,.... Anywhere $p1 shows up will be replaced with the text from the PARAMETER array. 
	TMP_FILENAME 
		This is the name of the directory the run will take place in. 
		******** WARNING *********
		Be sure to use all of the $p1, $p2,... in the file name so you don't overwrite data

	# echo ${TMP_FILENAME}
		this line is commented out, but it is used to test out file names. If you comment out the remaining lines of the for loop and run the script, the file names will list out to the console. If the look good, re-comment this line and uncomment the rest.
	Skip the next few lines to the python2.6 line.
		This line will genorate the pulse.inp file. By placing the $p1, $p2,... in the right locations you can change the various parameters in the input file. 
		******** WARNING *********
		the alph2 parameter defaults to 0.0d0 if it is not given. This is not the same as the default in the code.
		If you would like to see what is supported by the input_generator_2pulse.py code run 

		python input_generator_2pulse.py --help

		in the directory with the input_generator_2pulse.py file.
		******** WARNING *********
		Do not remove the > pulse.inp at the end of the line, since the code currently prints everything out to the console instead of a data file.
	the next line runs the code and puts the output into run.log
	Finally the cd .. moves you back out of the run directory and allows for the next run to happen.

titan.sh:
	The only major difference between titan.sh and stampede.sh is the 
	aprun -d16  
	in front of the running of the code.
	you will need to do the same update to the titan.sh as you did for the stampede.sh code. 

Let me know if you have any questions. There are also some tips in the comments of the .sh and .uni files. 