for i in $( ls 3pulse_rr2*.sh); do 
	echo $i;
	sbatch $i;
done
