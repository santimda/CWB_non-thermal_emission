set encoding iso_8859_1
set key top left
set xrange [-7:]
set yrange [26:]
set xtics 1
set xlabel 'log_{10} (E [eV])'
set ylabel 'log_{10} (E L(E) [erg/s])'

eV=1.602176e-12

set terminal postscript colour solid enhanced dl 1.5 lw 2.8 "Helvetica" 22
set output 'SED_syn.eps'

plot '../output/rad_syn1.dat' u (log10($1)):(log10($1*eV*$2)) t 'syn 1 unabs' w l lw 2 dt 1 lc 1,\
	'../output/rad_syn1.dat' u (log10($1)):(log10($1*eV*$3)) t 'syn 1 SSA' w l lw 2 dt 2 lc 1,\
	'../output/rad_syn2.dat' u (log10($1)):(log10($1*eV*$2)) t 'syn 2 unabs' w l lw 2 dt 1 lc 2,\
	'../output/rad_syn2.dat' u (log10($1)):(log10($1*eV*$3)) t 'syn 2 SSA' w l lw 2 dt 2 lc 2

!./fixbb SED_syn.eps
