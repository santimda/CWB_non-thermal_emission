set size ratio 2.
AU = 1.496e13
D = 80.0 # Stellar separation
set encoding iso_8859_1
set title 'L_{inj} in S1'
set xrange [0:1.5*D]
set yrange [0:2.2*D]
set cbrange [32:]
set xtics 20
set ylabel offset 1
set key off
set ylabel 'y [AU]'
set xlabel 'x [AU]'
set cblabel 'log(L_{inj} [erg s^{-1} AU^{-1}])'
set cbtics 0.5
set key off


# Stellar radii and positions (one is at 0,0; the other at D,0). 
# The 10*R is to make them bigger (i.e., not to scale)
R1=17.84*6.96e10/1.496e13
R2=13.06*6.96e10/1.496e13
# Plot the primary and secondary stars
set object 1 circle arc [0:90] at 0,0 size R1*10.0 front fc rgb "blue" fs solid 1.0
set label 1 "P" at D/10.,D/10. front textcolor rgb "blue" 
set object 2 circle arc [0:180] at D,0 size R2*10.0 front fc rgb "cyan" fs solid 1.0
set label 2 "S" at D*1.01,D/10. front textcolor rgb "cyan" 
#set arrow 1 from D/9, D/4 to D/25, D/25 head back filled linecolor rgb "blue" linewidth 2.1 
#set arrow 2 from D/1.2, D/4 to D/1.03, D/24. head back filled linecolor rgb "cyan" linewidth 2.1 

set terminal postscript color solid enhanced lw 2.5 "Helvetica" 20
set output 'Linj1_CD.eps'

plot '../output/Linj1.dat' u (($1)/AU):(($2)/AU):(log10(($3)/(($4)/AU))) w l lw 3 lc palette

!./fixbb Linj1_CD.eps
