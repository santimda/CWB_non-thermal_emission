set encoding iso_8859_1
set xrange [0.:1.]
unset key
set terminal postscript colour enhanced dl 2. lw 1.5 'Helvetica' 15
set output 'fluxes_phase.eps'
set multiplot layout 5,1  #rows, columns

set size ratio 0.4

#-------------------------------------------------------------------

# Enable the use of macros
set macros

# MACROS
# x- and y-tics for each row resp. column    
NOXTICS = "set xtics 0.2; set mxtics 2; set format x ''; unset xlabel"
XTICS = "set xtics ('0' 0.0, '' 0.1, '0.2' 0.2, '' 0.3, '0.4' 0.4, '' 0.5, '0.6' 0.6, '' 0.7, '0.8' 0.8, '' 0.9, '1' 1.0);\
      set mxtics 2; set xlabel '{/Symbol f}'"
NOYTICS = "set ylabel ' '"
YTICSL = "set ylabel 'F*10^{13} [erg s^{-1}]'; set ylabel offset 1.4,0"
YTICSR = "set ylabel 'F*10^{13} [erg s^{-1}]'; set ylabel offset 1.4,0"

set xtics nomirror

# Margins for each row 
VMARGIN1 = "set tmargin at screen 1.0; set bmargin at screen 0.80"
VMARGIN2 = "set tmargin at screen 0.8; set bmargin at screen 0.60"
VMARGIN3 = "set tmargin at screen 0.6; set bmargin at screen 0.40"
VMARGIN4 = "set tmargin at screen 0.4; set bmargin at screen 0.20"
VMARGIN5 = "set tmargin at screen 0.2; set bmargin at screen 0.00"

HMARGINL = "set lmargin at screen 0.04; set rmargin at screen 1.0"
HMARGINR = "set lmargin at screen 0.4; set rmargin at screen 1.0"

#-------------------------------------------------------------------

#set yrange [-0.95:1.4] 
#set ytics 0.5; set mytics 2

F0 = 1.0e-13

# First graphic (on top)
@VMARGIN1; @HMARGINL; @NOXTICS; @NOYTICS
LABEL = "> 0.1 TeV"
set obj 10 rect at graph 0.5,0.2 size char strlen(LABEL), char 2 
set obj 10 fillstyle solid border -1 front
set label 10 at graph 0.5,0.2 LABEL font "Helvetica,12" front center

plot '../output/fluxes_unabs.dat' u ($1):($7/F0) w l lw 2 lc 4 dt 2,\
	 '../output/fluxes_abs.dat' u ($1):($7/F0) w l lw 2 lc 1 dt 1 


# Second graphic
@VMARGIN2; @HMARGINL; @NOXTICS; @NOYTICS
LABEL = " 10-100 GeV "
set obj 10 rect at graph 0.5,0.2 size char strlen(LABEL), char 2 
set obj 10 fillstyle solid border -1 front
set label 10 at graph 0.5,0.2 LABEL font "Helvetica,12" front center

plot '../output/fluxes_unabs.dat' u ($1):($6/F0) w l lw 2 lc 4 dt 2,\
	 '../output/fluxes_abs.dat' u ($1):($6/F0) w l lw 2 lc 1 dt 1


# Third graphic
@VMARGIN3; @HMARGINL; @NOXTICS; @YTICSL
LABEL = "0.1-10 GeV"
set obj 10 rect at graph 0.5,0.2 size char strlen(LABEL), char 2 
set obj 10 fillstyle solid border -1 front
set label 10 at graph 0.5,0.2 LABEL font "Helvetica,12" front center

plot '../output/fluxes_unabs.dat' u ($1):($5/F0) w l lw 2 lc 4 dt 2,\
	 '../output/fluxes_abs.dat' u ($1):($5/F0) w l lw 2 lc 1 dt 1


# Fourth graphic
@VMARGIN4; @HMARGINL; @NOXTICS; @NOYTICS
LABEL = "0.1-100 MeV"
set obj 10 rect at graph 0.5,0.2 size char strlen(LABEL), char 2 
set obj 10 fillstyle solid border -1 front
set label 10 at graph 0.5,0.2 LABEL font "Helvetica,12" front center

plot '../output/fluxes_unabs.dat' u ($1):($4/F0) w l lw 2 lc 4 dt 2,\
	 '../output/fluxes_abs.dat' u ($1):($4/F0) w l lw 2 lc 1 dt 1


# Fifth graphic
@VMARGIN5; @HMARGINL; @XTICS; @NOYTICS
LABEL = "10-70 keV"
set obj 10 rect at graph 0.5,0.2 size char strlen(LABEL), char 2 
set obj 10 fillstyle solid border -1 front
set label 10 at graph 0.5,0.2 LABEL font "Helvetica,12" front center

plot '../output/fluxes_unabs.dat' u ($1):($3/F0) w l lw 2 lc 4 dt 2,\
	 '../output/fluxes_abs.dat' u ($1):($3/F0) w l lw 2 lc 1 dt 1


unset multiplot
!./fixbb fluxes_phase.eps
#!epstopdf Fluxes_ave.eps
