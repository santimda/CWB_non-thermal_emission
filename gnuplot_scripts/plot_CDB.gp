set encoding iso_8859_1
set title 'B evolution along the shock'
set key center right
set xrange [0:] # From 0 to ml/2
set yrange [:]
set xlabel 'l'
set ylabel 'B [mG]'

set terminal postscript color solid enhanced dl 6.5 lw 1.7 "Helvetica" 20
set output 'B_CD.eps'

plot '../output/thermo1.dat' u 1:($6*1.e3) t 'B_1' w l lw 2,\
	'../output/thermo2.dat' u 1:($6*1.e3) t 'B_2' w l lw 2

!./fixbb B_CD.eps

