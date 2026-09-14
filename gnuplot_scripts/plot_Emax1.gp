set encoding iso_8859_1
set title 'E_{max} in S1'
#set xtics 100
set ytics 1; set mytics 10 
set key top right
set ylabel 'log_{10}(E_{max} [eV])'
set xlabel 'position l'

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'Emax1.eps'

plot '../output/Emax_e1.dat' u (($1)):(log10(($2))) t 'e' w l lw 3,\
	'../output/Emax_p1.dat' u (($1)):(log10(($2))) t 'p' w l lw 3

!./fixbb Emax1.eps

