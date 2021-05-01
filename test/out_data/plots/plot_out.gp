set terminal pngcairo size 700,800
set output "out.png"

# set datafile separator ','
# set key autotitle columnhead

set multiplot

##### Constants ######
h = 6.626070040e-34
e = 1.6021766208e-19
c = 299792458

###### Plot parameters ######

set style line 1 lc rgb '#000000' lt 1 lw 2 pt 7 ps 1 dashtype 1 # black
set style line 2 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 1 # blue
set style line 3 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 1 # red
set style line 4 lc rgb '#006400' lt 1 lw 2 pt 6 ps 1 dashtype 1 # dark-green
set style line 5 lc rgb '#a52a2a' lt 1 lw 2 pt 8 ps 1 dashtype 1 # brown
set style line 6 lc rgb '#9400d3' lt 1 lw 2 pt 3 ps 1 dashtype 1 # purple
set style line 7 lc rgb '#008b8b' lt 1 lw 2 pt 5 ps 1 dashtype 1 # dark cyan
set style line 8 lc rgb '#ffd700' lt 1 lw 2 pt 7 ps 1 dashtype 1 # gold
set style line 9 lc rgb '#FF00ff' lt 1 lw 2 pt 7 ps 1 dashtype 1 # magenta
set style line 10 lc rgb '#000080' lt 1 lw 2 pt 7 ps 1 dashtype 1 # navy

set style line 21 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 2 # blue
set style line 31 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 2 # red

set style line 22 lc rgb '#0000ff' lt 1 lw 2 pt 2 ps 1 dashtype 4 # blue
set style line 32 lc rgb '#FF0000' lt 1 lw 2 pt 4 ps 1 dashtype 4 # red

set ylabel ""
set y2label ""

set mytics
set my2tics
set mxtics

set ytics nomirror # textcolor '#0000ff   '
set y2tics nomirror # textcolor '#FF0000'

set xrange[:]
set yrange[:]

set format y "% .3E"
set format y2 "% .3E"

#
##### Plot 1 #####
#

set origin 0,0.54
set size 1,0.46

set title "Electron density, U_{bias} = 0 meV, Gaussian barrier 100 meV"

set xlabel ""
set format x ""

set ylabel "n(x;k<0), n(x;k>0) [cm^{-3}]"
set y2label "n(x;k>0)-n(x;k<0) [cm^{-3}]"

set key over

plot '../../out' u 1:6 w l ls 2 title 'n(x;k<0)',\
    '' u 1:7 w l ls 3 title 'n(x;k>0)',\
    '' u 1:8 w l ls 4 axes x1y2 title 'n(x;k>0)-n(x;k<0)'

#
##### Plot 2 #####
#

set origin 0,0
set size 1,0.54

set title "Current density, U_{bias} = 0 meV, Gaussian barrier 20 meV"

set xlabel "x [nm]"
set format x "%.0f"
set ylabel "j(x;k<0), j(x;k>0) [Acm^{-2}]"
set y2label "j(x), j(x;k<0)+j(x;k>0) [Acm^{-2}]"

set key below

plot '../../out' u 1:9 w l ls 2 title 'j(x;k<0)',\
    '' u 1:10 w l ls 3 title 'j(x;k>0)',\
    '' u 1:11 w l ls 4 axes x1y2 title 'j(x;k<0)+j(x;k>0)',\
    '' u 1:3 w l ls 1 axes x1y2 title 'j(x)'
