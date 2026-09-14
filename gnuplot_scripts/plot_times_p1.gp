set encoding iso_8859_1
set size ratio 0.8
set title 'Proton cooling time (S1)'
ml = 300  # Number of points in the simulation
l = ml/10

set key top right
set xrange [9:]
set yrange [:16]
set ylabel 'log_{10}(t [s])'
set xlabel 'log_{10}(E_p [eV])'

set terminal postscript color solid enhanced dl 2.5 lw 2.0 "Helvetica" 20
set output 'tp1.eps'

# The cooling times depend on the position. Consider a fix position ($1 = X)
plot "< awk -v l=".l." '$1==l {print $2, $3}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'acc' w l lw 2.5 lc 8 dt 4,\
	"< awk -v l=".l." '$1==l {print $2, $4}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'pp' w l lw 2.5 lc 1,\
	"< awk -v l=".l." '$1==l {print $2, $5}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'adi' w l lw 2.5 lc 3,\
	"< awk -v l=".l." '$1==l {print $2, $6}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'conv' w l lw 2.5 lc rgb "coral" dt 6,\
	"< awk -v l=".l." '$1==l {print $2, $7}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'diff' w l lw 2.5 lc rgb "#32CD32" dt 6,\
	"< awk -v l=".l." '$1==l {print $2, $9}' ../output/CD_t_cool_p1.dat" u (log10($1)):(log10($2)) t 'cell' w l lc 7 dt 4

!./fixbb tp1.eps

