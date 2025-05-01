# ============================================================== #
# To call this file in gnuplot enter the following two commands:
# cd 'C:\Users\...'
# load 'C:\Users\...\plot_multi_beam_cond.gnuplot'
# ============================================================== #

### Set parameters to output .eps vector file with latex text:
set encoding utf8
set terminal epslatex size 7.0cm, 5.0cm
set output "multi_beam_cond.tex"
# set term png
# set output "multi_beam_cond.png"


### configure read in settings for .csv files
set datafile separator ","


### set axes
# set logscale x 10
set logscale y 10
set xrange [0.4 : 1.2]
set yrange [1E+06 : 1E+11]
# set grid


### label axes and graph
set xlabel '$o \ \text{(x-offset)}$' offset 0.0,0.5
set ylabel '$\kappa$' offset 1.5,0.0


### legend
set key at graph 0.71,0.975
# set key inside top center # position of the legend box
set key samplen 2.0 # make the sample lines in the legend box shorter
set key font ",10" # set the font size of the legend box
set key spacing 2.0 # set the spacing between the lines in the legend box
# set label 'a' at graph 0.87,0.3


### format axes
set format x "$%.1f$"
set format y "$10^{%L}$"
set xtics ( 0.4, 0.6, 0.8, 1.0, 1.2 ) offset 0.0,0.2
set ytics ( 1e6, 1e7, 1e8, 1e9, 1e10, 1e11 ) offset 0.53,0.0
#vertical spacing 1.0 samplen 1.1 width 1.5


### general plotting style
set style data linespoints


### set plot area
set lmargin at screen 0.16
set rmargin at screen 0.975
# set bmargin at screen 0.59
# set tmargin at screen 0.98


### plot
# use figure's color scheme: blue: #007AFF, red: #d20001, green: #7ec537, orange: #f59b23, dark grey: #636363, light grey: #9b9b9b
plot 'multi_beam_cond.csv' using 1:2 title 'Ghost: Off' pointtype 4 lw 3.5 ps 1.4 lt rgb "#007AFF" , \
     'multi_beam_cond.csv' using 1:3 title 'Ghost: On' pointtype 8 lw 3.5 ps 1.4 lt rgb "#d20001"


### write plots to files and close process
set output
# latex epslatex
set terminal pop
reset
cd 'C:\Users\...'


# ============================================================== #