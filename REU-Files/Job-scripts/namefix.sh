for i in $(ls); do
	if (echo $i | grep -q "_short_pulse_CEPA.uni"); then
		sed "s/_short_pulse/_short_pulse_CEPA/g"<$i | tee $i
	fi 
done