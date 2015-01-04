for i in $( ls ); 
do 
	if ( echo $i | grep -q "CEP-0.uni" ) ; then
		echo $i;
		sbatch $i;
	fi
done
echo done
