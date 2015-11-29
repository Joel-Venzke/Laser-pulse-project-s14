make 
# cp work-29-ne8-rr2-fixed run/
cd /Users/jvenzke/Laser_research/Data/H-skip_test
RUN_FILE=/Users/jvenzke/Laser_research/Laser-pulse-project-s14/REU-Files/Development/h_2h_skip/work-29-ne8-rr2-fixed


PARAMETER_1=( 1.0d0 1.1d0 1.25d0 1.5d0 2.0d0 0.90d0 0.95d0 )
PARAMETER_2=( 030.0d0 100.0d0 200.0d0 500.0d0 )

pwd
for p1 in ${PARAMETER_1[*]}; do
	for p2 in ${PARAMETER_2[*]}; do
		echo hfrac-$p1-xinside-$p2
		mkdir hfrac-$p1-xinside-$p2
		cd hfrac-$p1-xinside-$p2
		echo " &element    target = ' H'  /
 &discrete   nf = 6, nn(7)=4, nn(8)=4, nn(9)=4, nn(10)=4, ll(7)=0, ll(8)=1, ll(9)=2, ll(10)=3, 
             nn(11)=5, nn(12)=5, nn(13)=5, nn(14)=5, nn(15)=5, 
             ll(11)=0, ll(12)=1, ll(13)=2, ll(14)=3, ll(15)=4, nc = 4 /
 &control    key4 = 2, nprint = 25 /
 &numerics   dt=0.02d0, h=0.01d0, hfrac=$p1, xgrid=1.00d3, xinside=$p2, gbr=900.00d0, agbr=5.0d0/
 &energies   emin=0.002d0, de=0.002d0, nerg = 400 /" > tdse.inp
 		cat tdse.inp
 		cp /Users/jvenzke/Laser_research/Laser-pulse-project-s14/REU-Files/Development/h_2h_skip/run/pulse.inp .
 		$RUN_FILE> run.log
 		cd ..
	done
done