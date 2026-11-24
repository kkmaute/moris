# Run from a completed EXA run directory:
# gnuplot /path/to/projects/EXA/structure/dynamics/Toss_Ruler/plot_history.gp
# Optional: gnuplot -e "history='other.csv'; output_file='other.png'" plot_history.gp
if (!exists("history")) history = "Toss_Ruler_history.csv"
if (!exists("output_file")) output_file = "Toss_Ruler_history.png"
set datafile separator comma
set terminal pngcairo size 1400,1000 enhanced font "Sans,12"
set output output_file
set multiplot layout 2,2 title "Planar plane-stress toss ruler — DG(1) time integration"
set grid
set key left top
set xlabel "Time [s]"
set ylabel "Linear momentum [N s]"
plot history every ::1 using 1:4 with lines title "P_x", \
     history every ::1 using 1:5 with lines title "P_y"
set ylabel "Angular momentum about origin [N m s]"
plot history every ::1 using 1:6 with lines title "H_z"
set ylabel "Energy [J]"
plot history every ::1 using 1:7 with lines title "Kinetic", \
     history every ::1 using 1:8 with lines title "Strain", \
     history every ::1 using 1:9 with lines dashtype 2 title "Total"
set xlabel "Time [s]"
set ylabel "Right-midpoint displacement [m]"
plot history every ::1 using 1:11 with lines title "u_x", \
     history every ::1 using 1:12 with lines title "u_y"
unset multiplot
unset output
