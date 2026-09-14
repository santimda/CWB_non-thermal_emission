set encoding iso_8859_1
set key top left
set title 'thermal emission from WCR'
set xrange [0:6]
set yrange [24:34]
set xtics 3
set xlabel 'log (E [eV])'
set ylabel 'log (E*L(E) [erg/s])'

plot '../output/WCR_ff.dat' u (log10($1)):(log10($2)) t 'Thermal 1' w l,\
	'../output/WCR_ff.dat' u (log10($1)):(log10($3)) t 'Thermal 2' w l,\
	'../output/WCR_ff.dat' u (log10($1)):(log10($4)) t 'Thermal WCR' w l

set terminal postscript colour solid enhanced dl 2.5 lw 1.7 "Helvetica" 20
set output 'Lrad_WCRff.eps'
replot
