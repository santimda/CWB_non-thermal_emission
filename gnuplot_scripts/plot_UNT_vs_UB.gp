set encoding iso_8859_1
#set title '' U_{p,max} in FS'
#set xtics 100
#set ytics 1; set mytics 10 
set logscale
set key top right
set ylabel 'U_{NT}/U_{mag}'
set xlabel 'position l'

set terminal postscript color solid enhanced dl 6.5 lw 2.7 "Helvetica" 20
set output 'UNT_vs_UB.eps'

plot '../output/ENTp1.dat' u 1:4 t 'S1' w l lw 3,\
	 '../output/ENTp2.dat' u 1:4 t 'S2' w l lw 3 

#, '../output/ENTe_1.dat' u 1:3 w l lw 3  --> electrons are completely subdominant

!./fixbb UNT_vs_UB.eps

