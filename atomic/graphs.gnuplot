
set multiplot layout 1,2
set logscale x 
set xlabel 'r [au]'
set ylabel 'V_x(r)'
plot [1e-4:10] 'kli_pot.up' w l t 'V^{KLI}_x up', 'kli_pot.dw' w l t 'V^{KLI}_x down' 
set ylabel 'r V_x(r)'
plot [1e-4:15] 'kli_pot.up' u 1:($1*$2) w l t 'V^{KLI}_x up', 'kli_pot.dw'  u 1:($1*$2) w l t 'V^{KLI}_x down' 

#plot [1e-4:15] 'kli_pot.ud' u 1:($2-    $4) w l t 'V^{KLI}_x up - dw'
# unset multiplot
# set terminal postscript enhanced color
# set output 'kli.eps'
# replot
# unset output

pause -1