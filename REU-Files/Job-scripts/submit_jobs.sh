for i in $( ls 3pulse_rr2_ir*.sh); do 
	echo $i;
	sbatch $i;
done
