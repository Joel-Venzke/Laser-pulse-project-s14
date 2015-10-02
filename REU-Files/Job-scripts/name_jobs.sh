PAR=( 0.701d0 0.702d0 0.703d0 0.704d0 0.705d0 0.706d0 0.707d0 0.708d0 0.709d0 0.71d0 0.711d0 0.712d0 0.713d0 0.714d0 0.715d0 0.716d0 0.717d0 0.718d0 0.719d0 0.72d0 0.721d0 0.722d0 0.723d0 0.724d0 0.725d0 0.726d0 0.727d0 0.728d0 )

for i in ${PAR[*]}; do 
	# i2=$(echo $i | sed 's/3pulse_rr2/3pulse_rr2_ir/g')
	# cp 3pulse_rr2_ir_0.700d0.sh 3pulse_rr2_ir_$i.sh
	# subl 3pulse_rr2_ir_$i.sh
	# rm 3pulse_rr2_ir_$i.sh
done
