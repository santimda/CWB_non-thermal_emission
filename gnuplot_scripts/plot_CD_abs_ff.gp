set encoding iso_8859_1
set key top left
set title 'Free-free absorption'
ml = 300  # Number of points in the simulation
l = ml/2

set cbrange [-l:l]
set xrange [8:11]
set yrange [-0.05:1.05]
set xtics 1
set xlabel 'log_{10} ({/Symbol n} [Hz])'
set ylabel 'exp(-{/Symbol t})'
set cblabel 'position index (0 = apex)'

set terminal postscript colour solid enhanced lw 2. "Helvetica" 20
set output 'CD_absff.eps'

# The index was from 1 to ml, with 1 to ml/2 points for the top branch and from ml/2+1 to ml the bottom branch
plot "< awk -v l=".l." '$1<=l { print $1, $2, $3 }' ../output/CD_tauff.dat" u (log10($2)):(exp(-$3)):($1) t '' w l lc palette, \
	"< awk -v l=".l." '$1>l { print $1, $2, $3 }' ../output/CD_tauff.dat" u (log10($2)):(exp(-$3)):(l-$1) t '' w l lc palette

!./fixbb CD_absff.eps
