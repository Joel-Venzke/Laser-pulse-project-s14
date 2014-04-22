set term pslatex color 
set out 'approx_ionization_0375.tex'
set ylabel 'probability' offset -1.0,0.0
set xlabel 't [fs]' offset 0.0,-0.5
set format y '%3.1f'
#set mytics 10
set yr [0:1]
set ytics 0.2
set mytics 4
set xr [0.0:18]
#set xtics 0.05
set mxtics 5
#set format x '%4.2f'
set nolabel
set title '\textbf{2-36-2~~~$\mathbf{4 \times 10^{14}\,}$W/cm$^2$~~~S-S~~~10.2eV}'
set key center right spacing 1.3
plot datadir . '/S-s__2-36-2__1.0676d-1__0375/overlap.out' u ($1*.000120622):48 t '~~Ionization' w l lt 1 lc rgb "#009900" lw 2 , \
	datadir . '/S-s__2-36-2__1.0676d-1__0375/overlap.out' u ($1*.000120622):($37 + $38) t 'n $\le 2$' w l lt 1 lc 1 lw 2 , \
	datadir . '/S-s__2-36-2__1.0676d-1__0375/overlap.out' u ($1*.000120622):($47 - ($37 + $38)) t 'n $> 2$' w l lt 1 lc 3 lw 2 