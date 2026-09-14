set size ratio 2
AU = 1.496e13
D = 20 # 80.0 # Stellar separation
set xrange [0:]
set yrange [0:]
#set xtics 20
set key top left
set xlabel 'x [AU]'
set ylabel 'y [AU]'

# Stellar radii and positions (one is at 0,0; the other at D,0). 
# The 10*R is to make them bigger (i.e., not to scale)
R1 = 17.84*6.96e10/1.496e13
R2 = 13.06*6.96e10/1.496e13
# Plot the primary and secondary stars
set object 1 circle arc [0:90] at 0,0 size 10*R1 front fc rgb "blue" fs solid 1.0
set label 1 "Star 1" at D/25.,D/8. front 
set object 2 circle arc [0:180] at D,0 size 10*R2 front fc rgb "cyan" fs solid 1.0
set label 2 "Star 2" at 0.97*D,D/8. front 

set terminal postscript color solid enhanced dl 1.1 lw 2.0 "Helvetica" 20
set output 'CD.eps'

# Plot only the upper half of the CD
plot '../output/CD.dat' u ($1/AU):($2 >= 0 ? $2/AU : 1/0) t 'CD' w l lw 2 lc 3,\
	'../output/CD.dat' u ($1/AU):($2 >= 0 ? $2/AU : 1/0) t 'points' w p ps 0.5 pt 7 lc 7,\
 	'../output/CD_approx.dat' u ($1/AU):($2 >= 0 ? $2/AU : 1/0) t 'approx' w l lw 3 lc 2 dt 2
	
!./fixbb CD.eps
