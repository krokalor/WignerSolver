set terminal pngcairo size 800,500
set output "bc_comp.png"

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

# set title "TLE5-25"

set xlabel "p [a.u]"
set ylabel "f [a.u.]"

set mytics
set my2tics
set mxtics

# set ytics nomirror textcolor '#0000ff   '
# set y2tics nomirror textcolor '#FF0000'

set xrange[:]
set yrange[:4e-5]

set xtics -0.04,0.01,0.04

set key top right

# 0 mV
# set label at 9.15,3.52e-9 point pointtype 7 pointsize 1
# set label at 23.61,-7.897e-9 point pointtype 7 pointsize 1
# 10 mV
# set label at 5.63,1.00296e5 point pointtype 7 pointsize 1
# set label at 3.575,4.88959e5 point pointtype 7 pointsize 1
# 100 mV
# set label at 3.224,4.237e6 point pointtype 7 pointsize 1
# set label at 4.49,2.697e6 point pointtype 7 pointsize 1
# set label at 8.19,5.6134e6 point pointtype 7 pointsize 1

#
##### Plot #####
#

set origin 0,0
set size 0.5,1

plot 'BC_gauss.out' u 1:2 w l ls 1 axes x1y1 title 'Brak rozpraszań',\
    'BC_gauss_lorentz_1e-12.out' u 1:2 w l ls 2 axes x1y1 title 'Lorentz',\
    'BC_gauss_gauss_1e-12.out' u 1:2 w l ls 3 axes x1y1 title 'Gauss',\
    'BC_gauss_voigt_1e-12.out' u 1:2 w l ls 4 axes x1y1 title 'Voigt'

#
##### Inset #####
#

set origin 0.5,0
set size 0.5,1

set ylabel "f [a.u.]"

plot  'BC_gauss.out' u 1:2 w l ls 1 axes x1y1  title '{/Symbol t}_{/Symbol G} = 0 s',\
    'BC_gauss_lorentz_1e-11.out' u 1:2 w l ls 2 axes x1y1 title '{/Symbol t}_{/Symbol G} = 1E-11 s',\
    'BC_gauss_lorentz_1e-12.out' u 1:2 w l ls 3 axes x1y1  title '{/Symbol t}_{/Symbol G} = 1E-12 s',\
    'BC_gauss_lorentz_1e-13.out' u 1:2 w l ls 4 axes x1y1  title '{/Symbol t}_{/Symbol G} = 1E-13 s'
