# ============================================================== #
# To call this file in gnuplot enter the following two commands:
# Note: replace "..." with directory file is in (also at bottom of this script)
# cd 'C:\Users\...'
# load 'C:\Users\...\plot_multi_beam_H1s_error.gnuplot'
# ============================================================== #

### Set parameters to output .eps vector file with latex text:
set encoding utf8
set terminal epslatex size 6.2cm, 4.7cm
set output "multi_beam_L2_error.tex"
# set term png
# set output "multi_beam_L2_error.png"


### configure read in settings for .csv files
set datafile separator ","


### set axes
set logscale x 10
set logscale y 10
set xrange [1E+02 : 1E+06]
set yrange [1E-07 : 1E+00]
# set grid


### format axes
set format x "$10^{%L}$"
set format y "$10^{%L}$"
set xtics ( 1e2, 1e4, 1e6 ) offset 0.0,0.2
set ytics ( 1e-6, 1e-4, 1e-2, 1e0 ) offset 0.43,0.0
#vertical spacing 1.0 samplen 1.1 width 1.5


### label axes and graph
set xlabel '{\small Number of DOFs}' offset 0.0,0.6
set ylabel '{\footnotesize$\norm{\DiscDispl-\displ}{\mathrm{L2}}/\norm{\displ}{\mathrm{L2}}$}' offset 1.3,0.0


### legend
# set key at graph 0.69,0.975
set key inside top right # position of the legend box
set key samplen 15.0 # make the sample lines in the legend box shorter
set key font ",2" # set the font size of the legend box
set key spacing 2.9 # set the spacing between the lines in the legend box
# set label 'a' at graph 0.87,0.3


### general plotting style
set style data linespoints


### set plot area
set lmargin at screen 0.16
set rmargin at screen 0.975
# set bmargin at screen 0.59
set tmargin at screen 0.98


### plot
# use figure's color scheme: blue: #007AFF, red: #d20001, green: #7ec537, orange: #f59b23, dark grey: #636363, light grey: #9b9b9b
plot 'errors_non_enriched.csv' using 1:2 title 'Enrich.: Off ' pointtype 4 lw 3.5 ps 1.4 lt rgb "#007AFF" , \
     'errors_enriched.csv' using 1:2 title 'Enrich.: On ' pointtype 8 lw 3.5 ps 1.4 lt rgb "#d20001" , \
     'convergence.csv' using 1:2 title '$r=3$ ' with lines lw 2 dashtype 3 lt rgb "#636363"


### write plots to files and close process
set output
# latex epslatex
set terminal pop
reset
cd 'C:\Users\...'


# ============================================================== #