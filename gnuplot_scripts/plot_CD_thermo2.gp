set encoding iso_8859_1
set title 'Thermodynamical quantities along S2'
set key center right
set xrange [0:100] # From 0 to ml/2
set yrange [0:1.1] # The curves are normalized from 0 to 1
set xlabel 'l'
set ylabel ''

set terminal postscript color solid enhanced dl 3.5 lw 2.0 "Helvetica" 20
set output 'thermo_CD2.eps'

plot '../output/thermo_norm2.dat' u 1:2 t '{/Symbol r}/{/Symbol r}_{0}' w l lw 2,\
	'../output/thermo_norm2.dat' u 1:3 t 'P/P_0' w l lw 2,\
	'../output/thermo_norm2.dat' u 1:4 t 'v/v_0' w l lw 2,\
	'../output/thermo_norm2.dat' u 1:5 t 'T/T_0' w l lw 2,\
	'../output/thermo_norm2.dat' u 1:6 t 'B/B_{eq,0}' w l lw 2

!./fixbb thermo_CD2.eps
