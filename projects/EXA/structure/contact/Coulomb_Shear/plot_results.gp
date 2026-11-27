# Optional postprocessing only; no plotting package is needed by the EXA test.
# Usage: gnuplot -c plot_results.gp Coulomb_Shear_history.csv [output.png]
if (ARGC < 1) { print "Usage: gnuplot -c plot_results.gp history.csv [output.png]"; exit error }
set datafile separator ','
stats ARG1 every ::1 using 6 nooutput
force_scale = STATS_max
stats ARG1 every ::1 using 7 nooutput
slip_scale = STATS_max
set terminal pngcairo size 1600,640 enhanced font 'Sans,12'
set output (ARGC > 1 ? ARG2 : ARG1.'.png')
set multiplot layout 1,2 title 'Two-layer Coulomb shear cycle'
set grid
set xlabel 'D / D_{crit}'
set ylabel 'tau / (mu p)'
set key bottom right
plot ARG1 every ::1 using ($2/slip_scale):($6/force_scale) with lines dt 2 lc rgb 'black' title 'Analytic endpoints', \
     ARG1 every ::1 using ($2/slip_scale):($3/force_scale) with linespoints pt 7 ps 0.7 title 'MORIS Seitz'
set ylabel 's / D_{crit}'
set key top left
plot ARG1 every ::1 using ($2/slip_scale):($7/slip_scale) with lines dt 2 lc rgb 'black' title 'Analytic endpoints', \
     ARG1 every ::1 using ($2/slip_scale):($5/slip_scale) with linespoints pt 7 ps 0.7 title 'MORIS Seitz'
unset multiplot
unset output
