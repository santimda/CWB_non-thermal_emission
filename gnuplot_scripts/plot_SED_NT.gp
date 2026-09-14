set encoding iso_8859_1
set key top left
set title '' #SED for electrons and protons
set xrange [-7:14]
set yrange [26:34]
set xtics 3; set mxtics 3
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'

eV=1.602176e-12

set terminal postscript colour solid enhanced dl 1.5 lw 1.7 "Helvetica" 20
set output 'SED_NT.eps'

plot '../output/rad_ic_WCR.dat' u (log10($1)):(log10($2*$1*eV)) t 'unabs IC' w l lw 2 dt 2 lc 1,\
	'../output/rad_ic_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'abs IC' w l lw 2 lc 1,\
	'../output/rad_br_WCR.dat' u (log10($1)):(log10($2*$1*eV)) t 'unabs Br' w l lw 2 dt 2 lc 2,\
	'../output/rad_br_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'abs Br' w l lw 2 lc 2,\
	'../output/rad_syn_WCR.dat' u (log10($1)):(log10($2*$1*eV)) t 'Sync + R-T' w l lw 2 dt 2 lc 4,\
	'../output/rad_syn_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'Sync + R-T + FFA' w l lw 2 lc 4,\
	'../output/rad_pp_WCR.dat' u (log10($1)):(log10($2*$1*eV)) t 'unabs pp' w l lw 2 dt 2 lc 5,\
	'../output/rad_pp_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'abs pp' w l lw 2 lc 5

!./fixbb SED_NT.eps
