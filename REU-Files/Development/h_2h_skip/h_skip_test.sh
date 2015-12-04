#!/bin/bash
#!/bin/bash
#PBS -N supermic
#PBS -j oe
#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=20

# run the shell script, sit back, have a sip of coffee, and read this https://xkcd.com/676/
export OMP_NUM_THREADS=20
export KMP_STACKSIZE=16m

RUN_DIR=/work/jvenzke/H-skip_test/
mkdir -p $RUN_DIR
cd $RUN_DIR
RUN_FILE=/home/jvenzke/LaserPulseShapeStudies/REU-Files/Development/h_2h_skip/work-29-ne8-rr2-fixed
PULSE_FILE=/home/jvenzke/LaserPulseShapeStudies/REU-Files/Development/h_2h_skip/pulse.inp
NE_FILE=/home/jvenzke/LaserPulseShapeStudies/REU-Files/Development/Input-files/ne.wfn

PARAMETER_1=( 1.0d0 1.1d0 1.25d0 1.5d0 2.0d0 0.90d0 0.95d0 )
PARAMETER_2=( 030.0d0 100.0d0 200.0d0 500.0d0 )

for p1 in ${PARAMETER_1[*]}; do
	for p2 in ${PARAMETER_2[*]}; do
 		mkdir hfrac-$p1-xinside-$p2
		cd hfrac-$p1-xinside-$p2
		echo " &element    target = 'Ne'  /
 &discrete   nf=26, nnauto = 2, llauto=1, mfixed=0, nc=4, nn(1)=1, ll(1)=0,
             nn(2)=2, ll(2)=0, nn(3)=3, ll(3)=0, nn(4)=4, ll(4)=0, nn(5)=5, ll(5)=0,
             nn(6)=6, ll(6)=0, nn(7)=7, ll(7)=0, nn(8)=8, ll(8)=0, 
             nn(9)=2, ll(9)=1, nn(10)=3, ll(10)=1, nn(11)=4, ll(11)=1, nn(12)=5, ll(12)=1,
             nn(13)=6, ll(13)=1, nn(14)=7, ll(14)=1, nn(15)=8, ll(15)=1, 
             nn(16)=3, ll(16)=2, nn(17)=4, ll(17)=2, nn(18)=5, ll(18)=2,
             nn(19)=6, ll(19)=2, nn(20)=7, ll(20)=2, nn(21)=8, ll(21)=2,
             nn(22)=4, ll(22)=3, nn(23)=5, ll(23)=3,
             nn(24)=6, ll(24)=3, nn(25)=7, ll(25)=3, nn(26)=8, ll(26)=3 /
 &control    key4 = 2, nprint = 500, irestart = 0/
 &numerics   dt = 0.0055d0, h = 0.02d0, hfrac=$p1, xgrid=3.20d3, xinside=$p2, gbr = 3000.00d0, agbr = 5.0d0/
 &energies   emin = 0.45d0, de = 0.0001d0, nerg = 3750 /">tdse.inp
 		cp $NE_FILE .
		cp $PULSE_FILE .
 		${RUN_FILE} > run.log
 		cd ..
	done
done
