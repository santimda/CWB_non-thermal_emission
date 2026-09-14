set encoding iso_8859_1
set title 'Synthetic map at 2.3 GHz (abs)'
set key off
set xlabel 'x [mas]'
set ylabel 'y [mas]'
set cblabel 'Flux [mJy]'

D = 47.0
a = 2.*D
set yrange [-0.5*a:1.5*a]
set xrange [-1*a:1*a]

set cbrange [-3.0:-1.0]
set size ratio 2./2.

set label 1 "" at +0.002*D,0 front point pointtype 1 lc 3 lw 2 ps 16
set label 2 "1" at 0.25*D,-a/10. font "Helvetica,23" front textcolor lt 3
set label 3 "" at -0.002*D,D front point pointtype 1 lc 1 lw 2 ps 16
set label 4 "2" at 0.25*D,1.25*D font "Helvetica,23" front textcolor lt 1

set pm3d explicit
set view map

splot '../output/WCR_radiomap_2GHz.dat' u 2:(-$1):(log10($4)) w points palette
set terminal postscript colour enhanced font "Helvetica,20"
set output 'radiomap_2GHz_abs.eps'
replot
!./fixbb radiomap_2GHz_abs.eps
