set encoding iso_8859_1
set key top left
set title 'SSA absorption in S1'
ml = 300  # Number of points in the simulation
l = ml/2

set cbrange [0:l]
set xrange [7.7:10]
set yrange [-0.01:1.01]
set xtics 1
set xlabel 'log_{10} ({/Symbol n} [Hz])'
set ylabel '(1 - exp(-{/Symbol t}))/{/Symbol t}'
set cblabel 'index l (distance from apex)'

set terminal postscript colour solid enhanced lw 2. "Helvetica" 20
set output 'opac_SSA1.eps'

#Upper half
plot "< awk -v l=".l." '$1<=l { print $1, $2, $3 }' ../output/opac_SSA1.dat" u (log10($2)):3:1 t '' w l lc palette

!./fixbb opac_SSA1.eps
