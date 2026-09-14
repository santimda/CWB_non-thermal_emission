eV=1.602176e-12
set encoding iso_8859_1
set title ''
#'Electron distribution in SW1'
set xrange [6:13.5]
set yrange [30:]
set xlabel 'log_{10}(E_e [eV])'; set ylabel offset 0.0,-0.1
set ylabel 'log_{10}({E_e}^2 N_e(E_e) [erg])'; set ylabel offset 0.1,0.0
set ytics 1 ; set mytics 2
set mxtics 2
set palette maxcolors 5 # 5 = N_bins
set cblabel '' #'{/Symbol q}/{/Symbol p}'; set cblabel offset -0.5,-0.8
set cbtics ("I1" 1.4, "I2" 2.2,"I3" 3.,"I4" 3.8,"I5" 4.6)
set palette defined ( 0 "#FF8C00", 1 "#8B0000" ) # From dark-orange to dark-red
set datafile missing

set terminal postscript color solid enhanced dl 1.5 lw 2.5 "Helvetica" 21
set output 'dist_e_wcr1_bin.eps'

plot '../output/dist_e_wcr1_bin.dat' u (log10($2)):(log10((($2)*eV)**2*($3))):1 w l lw 3.0 dt 1 lc palette t '',\
     '../output/dist_e1_total.dat' u (log10($1)):(log10((($1)*eV)**2*($2))) w l lw 3.0 dt 2 lc 0 t 'Total'	 

!./fixbb dist_e_wcr1_bin.eps
