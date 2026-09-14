set encoding iso_8859_1
set title '' #'Electron distribution in S1'
ml = 300  # Number of points in the simulation
set cbrange [0:ml/2]
set xrange [:]
set yrange [26:]
set cbrange [:] # from 1 to ml/2 
set xlabel 'log_{10}(E_e [eV])'; set ylabel offset 0.0,-0.1
set ylabel 'log_{10}({E_e}^2 N_e(E_e) [erg])'; set ylabel offset 0.1,0.0
set cblabel 'position (l)'; set cblabel offset -0.2,-1.1
set key off
set datafile missing

eV=1.602176e-12

set terminal postscript color solid enhanced dl 2. lw 2. "Helvetica" 22
set output 'dist_e1_CD.eps'

# 1 < $1 <= ml/2
plot '../output/dist_e1.dat' \
     u (($1 > 1 && $1 <= ml/2) ? log10($2) : 1/0):(($1 > 1 && $1 <= ml/2) ? log10(($2*eV)**2.*$3) : 1/0):1 \
     w l lw 1.0 lc palette, \
     '../output/dist_e1_total.dat' \
     u (log10($1)):(log10(($1*eV)**2*$2)) w l lw 3.0 dt 2

!./fixbb dist_e1_CD.eps
