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
a = 4.0*pi*dist**2
h = 6.6260755e-27
eV = 1.602176e-12

# Plot the broadband SED. Fit each NT emission component to a function and get the total SED as the sum.

x0=8.0
f(x)=a0+a1*(x-x0)+a2*(x-x0)**2+a3*(x-x0)**3+a4*(x-x0)**4+a5*(x-x0)**5+a6*(x-x0)**6+a7*(x-x0)**7+a8*(x-x0)**8
fit f(x) "< awk '($1 >8.e0 && $1<5.e12) { print $1, $2, $3 }' ../output/rad_ic_WCR.dat" u (log10($1)):(log10($3*$1*eV)) via a0,a1,a2,a3,a4,a5,a6,a7,a8
f1(x)=x>13?0 : x<0.5?0 : f(x)

x1=-1.0
g(x)=b0+b1*(x-x1)+b2*(x-x1)**2+b3*(x-x1)**3+b4*(x-x1)**4+b5*(x-x1)**5+b6*(x-x1)**6+b7*(x-x1)**7
fit g(x) "< awk '($1 >1.e-6 && $1<1.e4) { print $1, $2, $3 }' ../output/rad_syn_WCR.dat" u (log10($1)):(log10($3*$1*eV)) via b0,b1,b2,b3,b4,b5,b6,b7
g1(x)=x>6?0 : x<-6?0 : g(x)

x2=7.0
j(x)=c0+c1*(x-x2)+c2*(x-x2)**2+c3*(x-x2)**3+c4*(x-x2)**4+c5*(x-x2)**5+c6*(x-x2)**6+c7*(x-x2)**7
fit j(x) "< awk '($1 >1.e5 && $1<1.e12) { print $1, $2, $3 }' ../output/rad_br_WCR.dat" u (log10($1)):(log10($3*$1*eV)) via c0,c1,c2,c3,c4,c5,c6,c7
j1(x)=x>12?0 : x<5?0 : j(x)

x3=11.0
k(x)=d0+d1*(x-x3)+d2*(x-x3)**2+d3*(x-x3)**3+d4*(x-x3)**4+d5*(x-x3)**5+d6*(x-x3)**6+d7*(x-x3)**7
fit k(x) "< awk '($1 >1.e8 && $1<1.e14) { print $1, $2, $3 }' ../output/rad_pp_WCR.dat" u (log10($1)):(log10($3*$1*eV)) via d0,d1,d2,d3,d4,d5,d6,d7
k1(x)=x>14?0 : x<8?0 : k(x)


set terminal postscript colour solid enhanced dl 6.5 lw 2.8 "Helvetica" 22
set output 'SED_NT_sum.eps'

# Plot each component, some data points and the total emission from the model
plot '../output/rad_ic_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t '' w l lw 1,\
     '../output/rad_syn_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t '' w l lw 1,\
     '../output/rad_pp_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t '' w l lw 1,\
     '../output/rad_br_WCR.dat' u (log10($1)):(log10($3*$1*eV)) t '' w l lw 1,\
     '../obs_data/data_radio2008y9.dat' u (log10(($1)*1.0e9*h/eV)):(log10(($1)*1.0e9*a*($2)*1.0e-26)):(log10(1+($3)/($1))):(log10(1.+($4)/($2))) w xyerrorbars pt 7 t '',\
     '../obs_data/sens_CTA_100h.dat' u (log10(1.e6*($1))):(log10(($2)*1.6*a)) t 'CTA' w l lt 2 lc rgb "black" lw 5,\
     '../obs_data/sens_fermi_0,90.dat' u (log10($1)):(log10(($2)*a)) t 'Fermi' w l lt 2 lc rgb '#77252525' lw 5,\
      log10(10**f1(x)+10**g1(x)+10**j1(x)+10**k1(x)) t 'model' lw 3


!./fixbb SED_NT_sum.eps
