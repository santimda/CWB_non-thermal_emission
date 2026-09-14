set encoding iso_8859_1
set key top left
set xrange [-6:14]
set yrange [26:35]
set xtics 3
set mxtics 3
set ytics 2
set mytics 2
set xlabel 'log_{10} (E [eV])'
set ylabel 'log_{10} (E L(E) [erg/s])'
set pointsize 1.5

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 est� en GHz y Sv en $2 en mJy
dist = 2.4*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
a = 4.0*pi*dist**2

# Plot each NT emission component, some data points and sensitivity curves
plot '../output/rad_ic_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'IC' w l lw 2,\
	'../output/rad_syn_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'Sync' w l lw 2,\
	'../output/rad_pp_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'p-p' w l lw 2,\
	'../output/rad_br_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'Bre' w l lw 2,\
	'../obs_data/data_radio2008y9.dat' u (log10(($1)*10**9*h/eV)):(log10(($1)*10**9*a*($2)*10**(-26))):(log10(1+($3)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 t '',\
	'../obs_data/sens_CTA_100h.dat' u (log10(1.e9*($1))):(log10(($2)*1.6*a)) t 'CTA' w l lt 2 lc rgb "black" lw 5,\
	'../obs_data/sens_fermi_0,90.dat' u (log10($1)):(log10(($2)*a)) t 'Fermi' w l lt 2 lc rgb '#77252525' lw 5

set terminal postscript colour solid enhanced dl 6.5 lw 2.8 "Helvetica" 22
set output 'modelNT.eps'
replot
!./fixbb modelNT.eps
