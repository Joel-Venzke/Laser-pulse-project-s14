make 
cp work-29-ne8-rr2-fixed run/
cd run
./work-29-ne8-rr2-fixed

files=( hyd1.out hyd3.out hydtotion.out )
stub="test/"
# stub="short/"
for i in ${files[*]}; do
	echo =========================================================================================
	echo =========================================================================================
	echo =========================================================================================
	echo $i
	echo =========================================================================================
	echo =========================================================================================
	echo =========================================================================================
	diff $i $stub$i
	# read -p "press enter"
done