for i in $( ls ); 
do 
	if ( echo $i | grep -q ".uni" ) ; then
		echo $i;
		sbatch $i;
	fi
done
