for i in $( ls ); 
do 
	if ( echo $i | grep -q "_pulse_CEPA.uni" ) ; then
		echo $i;
		sbatch $i;
	fi
done
