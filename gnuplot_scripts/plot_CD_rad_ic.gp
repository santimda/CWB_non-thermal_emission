set encoding iso_8859_1
set key top left
set xrange [1:]
set yrange [:]
set xtics 1
set xlabel 'log_{10} (E [eV])'
set ylabel 'log_{10} (E L(E) [erg/s])'

set terminal postscript colour solid enhanced dl 1.5 lw 2.8 "Helvetica" 22
set output 'SED_ic.eps'

plot '../output/rad_ic1.dat' u (log10($1)):(log10($2)) t 'IC 1 unabs' w l lw 2 dt 2 lc 1,\
	'../output/rad_ic1.dat' u (log10($1)):(log10($3)) t 'IC 1 abs' w l lw 2 dt 1 lc 1,\
	'../output/rad_ic2.dat' u (log10($1)):(log10($2)) t 'IC 2 unabs' w l lw 2 dt 2 lc 2,\
	'../output/rad_ic2.dat' u (log10($1)):(log10($3)) t 'IC 2 abs' w l lw 2 dt 1 lc 2

!./fixbb SED_ic.eps
