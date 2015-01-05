for i in $( ls ); 
do 
	if ( echo $i | grep -q "short_pulse.uni" ) ; then
		echo $i;
		sbatch $i;
	fi
done
echo done
