set encoding iso_8859_1
set key top left
set xrange [-7:15]
set yrange [25:33]
set xtics 3
set xlabel 'log (E [eV])'
set ylabel 'log (E*L(E) [erg/s])'
eV = 1.602176e-12

plot '../output/rad_syn1.dat' u (log10($1)):(log10($2*$1*eV)) t 'syn 1' w l lw 2, '../output/rad_syn2.dat' u (log10($1)):(log10($2*$1*eV)) t 'syn 2' w l lw 2,'../output/rad_ic1.dat' u (log10($1)):(log10($3)) t 'IC 1' w l lw 2,'../output/rad_ic2.dat' u (log10($1)):(log10($3)) t 'IC 2' w l lw 2, '../output/rad_br1.dat' u (log10($1)):(log10($3)) t 'Brems 1' w l lw 2, '../output/rad_br2.dat' u (log10($1)):(log10($3)) t 'Brems 2' w l lw 2, '../output/rad_pp1.dat' u (log10($1)):(log10($3)) t 'p-p 1' w l lw 2, '../output/rad_pp2.dat' u (log10($1)):(log10($3)) t 'p-p 2' w l lw 2

#'./rad_pp1.dat' u (log10($1)):(log10($2)) t 'p-p 1 (unabs)' w l lw 2, './rad_pp1.#dat' u (log10($1)):(log10($3)) t 'p-p 1 (abs)' w l lw 2
set terminal postscript colour solid enhanced dl 6.5 lw 1.7 "Helvetica" 20
set output 'Lrad.eps'
replot
