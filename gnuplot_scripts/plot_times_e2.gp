set encoding iso_8859_1
set size ratio 0.8
set title 'Electron cooling time (S2)'
ml = 300  # Number of points in the simulation
l = ml/10

set key box opaque height 1.01 width 0.50 
set key samplen 2.5 spacing 1.01 font ",20"
set key bottom left
set xrange [6.0:13]
set yrange [-1:12]
set ylabel 'log_{10}(t [s])'; set ylabel offset 0.8,0
set xlabel 'log_{10}(E_e [eV])'; set xlabel offset 0,0.4
set ytics 2; set mytics 2

set terminal postscript color solid enhanced dl 1.5 lw 2. "Helvetica" 20
set output 'te2.eps'

# The cooling times depend on the position.
# Consider a fixed position at l=ml//10 ($1 = l)
plot "< awk -v l=".l." '$1==l {print $2, $9 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'acc' w l lw 2.5 lc 8 dt 4,\
    "< awk -v l=".l." '$1==l {print $2, $5 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'IC' w l lw 2.5 lc 6,\
    "< awk -v l=".l." '$1==l {print $2, $6 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'sy' w l lw 2.5 lc rgb "#DC143C",\
    "< awk -v l=".l." '$1==l {print $2, $7 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'br' w l lw 2.5 lc 5,\
    "< awk -v l=".l." '$1==l {print $2, $8 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'diff' w l lw 2.5 lc rgb "#32CD32" dt 6,\
    "< awk -v l=".l." '$1==l {print $2, $11 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'conv' w l lw 2.5 lc rgb "coral" dt 6,\
    "< awk -v l=".l." '$1==l {print $2, $12 }' ../output/CD_t_cool_e2.dat" u (log10($1)):(log10($2)) t 'cell' w l lc 7 dt 4

!./fixbb te2.eps