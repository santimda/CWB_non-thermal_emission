set view map

frac = 0.8 # ~60 mJy / 80 mJy ?

ang_deg = 185.#275.0
ang = ang_deg*pi/180.  # rotation angle to match observed map (Marcote2021)

dec_deg = -51.7 
cos_dec = cos(dec_deg * pi/180.)

set contour base
set isosamples 250
set cntrparam bspline
#set cntrparam level incremental -2.,0.25,1.0
set cntrparam levels discrete 0.114, 0.16122, 0.228, 0.32244, 0.456, 0.64488, 0.912, 1.28976, 1.824, 2.57952, 3.648, 5.15905, 7.296, 10.3181, 14.592
unset surface
set table 'cont.dat'
splot '../output/WCR_conv_2Gabs.dat' u ($1*cos(ang)+$2*sin(ang)):($1*sin(ang)-$2*cos(ang)):($3*frac)
unset table

reset
set encoding iso_8859_1
#set title 'Synthetic map at 2.3 GHz'
set title ''

D = 47.0
a = 2.*D
set xlabel 'x [mas]'
set ylabel 'y [mas]'
#set cblabel 'log_{10}(Flux [mJy beam^{-1}])'
set cblabel 'Flux [mJy beam^{-1}]'
set yrange [-60:60]
#set xrange [-100:50]
set xrange [-90:30]
set size ratio 120/120 #(120*cos_dec)
set cbrange [0.0:10]
set logscale cb


# Position of the stars: WN at (0,0), WC at (D,0) --> rotate the second one
Dx = D*cos(ang)
Dy = D*sin(ang)

set pm3d interpolate 0,0
set view map
#set palette gray negative
set palette defined (0 "white", 1 "grey50")
unset key
set label 1 "" at 0,0 front point pointtype 1 lc 3 lw 2 ps 8
set label 2 "WN" at 0.2*D,-a/10. font "Helvetica,23" front textcolor lt 3

set label 3 "" at Dx,Dy front point pointtype 1 lc 3 lw 2 ps 8
set label 4 "WC" at 1.4*Dx,-0.7*Dy font "Helvetica,23" front textcolor lt 3

siz_x = 11.3
siz_y = 5.6
set object 5 ellipse center graph 0.11,0.11 size siz_y,siz_x angle ang fc rgb "grey30" fs solid 1.0 noborder front

splot '../output/WCR_conv_2Gabs.dat' u ($1*cos(ang)+$2*sin(ang)):($1*sin(ang)-$2*cos(ang)):($3*frac) w pm3d,\
      'cont.dat' w l lt -1 lw 1.5

set terminal postscript colour enhanced dl 2.5 lw 1.7 "Helvetica" 22
set output 'convmap_2GHz.eps'
replot
!./fixbb convmap_2GHz.eps
#set output 'convmap_2GHz_i10_nohard_1D.eps'
#replot
#!./fixbb convmap_2GHz_i10_nohard_1D.eps
