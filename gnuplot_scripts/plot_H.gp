set encoding iso_8859_1
set size ratio 0.85
set title 'Shock width'
#set xtics 100
#set xrange [0.:2.4]
#set yrange [-0.1:0.5]
#set logscale y
#set ytics 1; set mytics 10
#set xtics 0.5; set mxtics 5; set xtics nomirror font ",16" 
set ylabel 'H/R'
set xlabel 'l'; set xlabel offset 0,0.25 font ",18"
#set label "NAME" at graph 0.05,0.95 font "Helvetica,18"

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'H.eps'

plot '../output/H1.dat' u 1:4 t 'S1' w l lw 2.4 dt 1 lc 1,\
     '../output/H2.dat' u 1:4 t 'S2' w l lw 2.4 dt 1 lc 2

!./fixbb H.eps