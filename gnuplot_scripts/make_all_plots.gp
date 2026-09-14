#! /bin/bash

gnuplot <<-EOF
load 'plot_abs_SSA1.gp'
reset
load 'plot_CD.gp'
reset
load 'plot_CD_abs_ff.gp'
reset
load 'plot_CD_dist_e1.gp'
reset
load 'plot_CD_dist_e2.gp'
reset
load 'plot_CD_dist_p1.gp'
reset
load 'plot_CD_dist_p2.gp'
reset
load 'plot_CD_Linj1.gp'
reset
load 'plot_CD_Linj2.gp'
reset
load 'plot_CD_rad.gp'
reset
load 'plot_CD_rad_br.gp'
reset
load 'plot_CD_rad_ic.gp'
reset
load 'plot_CD_rad_pp.gp'
reset
load 'plot_CD_rad_syn.gp'
reset
load 'plot_CD_thermo1.gp'
reset
load 'plot_CD_thermo2.gp'
reset
load 'plot_times_e1.gp'
reset
load 'plot_times_e2.gp'
reset
load 'plot_times_p1.gp'
reset
load 'plot_times_p2.gp'
reset
load 'plot_CDB.gp'
reset
load 'plot_dist_e_wcr1_bin.gp'
reset
load 'plot_dist_e_wcr2_bin.gp'
reset
load 'plot_dist_p_wcr1_bin.gp'
reset
load 'plot_dist_p_wcr2_bin.gp'
reset
load 'plot_Emax1.gp'
reset
load 'plot_Emax2.gp'
reset
load 'plot_H.gp'
reset
load 'plot_SED_NT.gp'
reset
load 'plot_UNT_vs_UB.gp'
reset
load 'plot_WCR_ff.gp'
reset
load 'plot_WCR_rad.gp'
quit
EOF
