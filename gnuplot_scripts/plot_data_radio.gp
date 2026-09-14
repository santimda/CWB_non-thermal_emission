set encoding iso_8859_1
set key bottom right
set xrange [-6.0:-3.5]
#set yrange [28.:30.5]
set xtics 1; set mxtics 10
set ytics 1; set mytics 10
set xlabel 'log_{10} ({/Symbol e} [eV])'
set ylabel 'log_{10} ({/Symbol e} L({/Symbol e}) [erg s^{-1}])'
set pointsize 2
set size square

# 1 Jy = 10^-23 erg/cm2/s/Hz, nu en $1 est� en GHz y Sv en $2 en mJy
dist = 2.4*3.0856e21
h = 6.6260755e-27
eV = 1.602176e-12
f = h/eV
a = 4.0*pi*dist**2

# Fit the thermal and non-thermal components and them add them in a new function
x0=-5.0
fs(x)=a0+a1*(x-x0)+a2*(x-x0)**2+a3*(x-x0)**3+a4*(x-x0)**4+a5*(x-x0)**5+a6*(x-x0)**6
fit fs(x) "< awk '($1 >8.e-7 && $1<1.e-3) { print $1, $2, $3 }' ../output/rad_syn_WCR.dat" u (log10($1)):(log10($3*$1*eV)) via a0,a1,a2,a3,a4,a5,a6

fw(x) = aw*x+bw
fit fw(x) '../output/windemiff.dat' u (log10((($1)*h/eV))):(log10(($1)*a*($2)*10**(-26))) via aw,bw

f(x) = log10( 10.**(fs(x)) + 10.**(fw(x)) )

set terminal postscript colour solid enhanced dl 1.5 lw 2.8 "Helvetica" 22
set output 'dataradio_vs_model.eps'

# Plot the data points, each curve and the sum
plot '../output/rad_syn_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t 'Sync' w l lw 2,\
    '../output/windemiff.dat' u (log10($1*f)):(log10(($1)*a*($2)*1.0e-26)) w l lw 2 lc 7 t 'wind f-f',\
    '../output/windemiff1.dat' u (log10($1*f)):(log10(($1)*a*($2)*1.0e-26)) w l dt 2 lw 1 t 'wind1 f-f',\
    '../output/windemiff2.dat' u (log10($1*f)):(log10(($1)*a*($2)*1.0e-26)) w l dt 4 lw 1 t 'wind2 f-f',\
    f(x) t 'total' lw 3 lc -1 dt 3

!./fixbb SED_radio.eps
