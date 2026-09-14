set encoding iso_8859_1
set key box opaque top right
set xrange [8.0:11]  
set yrange [0:3]   # [0:3.5]
set xtics 1; set mxtics 10
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol n} [Hz])'
set ylabel 'log_{10} (S_{/Symbol n} [mJy])'
set pointsize 2
#set size square  # Good for two-column texts

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 est� en GHz y Sv en $2 en mJy
mJy = 1.0e-26
dist = 2.4*3.0856e21  	# Distance to the source = 2.4 kpc
h = 6.6260755e-27 		#erg.seg
eV = 1.602176e-12 		# 1 eV = 1.602176e-12 erg
a = 4.0*pi*dist**2

# Fit the thermal and non-thermal components and them add them in a new function
x0=-5.5
fsy(x)=a0+a1*(x-x0)+a2*(x-x0)**2+a3*(x-x0)**3+a4*(x-x0)**4+a5*(x-x0)**5+a6*(x-x0)**6+a7*(x-x0)**7
fit fsy(x) "< awk '($1 >3.e-7 && $1<8.e-4) { print $1, $2, $3 }' ../output/rad_syn_WCR.dat" u (log10($1*(eV/h))):(log10($3*(h/a/mJy))) via a0,a1,a2,a3,a4,a5,a6,a7

fw(x) = aw*(x-x0)+bw
fit fw(x) '../output/windemiff.dat' u (log10($1)):(log10($2)) via aw,bw

f(x) = log10( 10.**(fsy(x)) + 10.**(fw(x)) )

# Plot the data points, each curve and the sum
plot '../python_scripts/data_radio_apep.dat' u (log10($1*1e9)):(log10($2)):(log10(1+$3/$1)):(log10(1+$4/$2)) w xyerrorbars pt 7 t 'Data',\
    '../output/rad_syn_WCR.dat' u (log10($1*(eV/h))):(log10($3*(h/a/mJy))) t 'sync' w l lw 1,\
    '../output/windemiff1.dat' u (log10($1)):(log10($2))  w l dt 2 lw 1 t 'wind1 f-f',\
    '../output/windemiff2.dat' u (log10($1)):(log10($2)) w l dt 3 lw 1 t 'wind2 f-f',\
    f(x) t 'total' lw 2


set terminal postscript colour solid enhanced dl 6.5 lw 2.8 "Helvetica" 22
set output 'dataradio_vs_model.eps'
replot
!./fixbb dataradio_vs_model.eps
