for i in $( ls 3pulse_rr2*.sh); do 
	echo $i;
	i2=$(echo $i | sed 's/3pulse_rr2/3pulse_rr2_ir/g')
	cp $i $i2
done
