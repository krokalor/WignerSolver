set terminal pngcairo enhanced size 500,500
set output 'ivChar.png'

set multiplot

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

#
# Plot
#

set ylabel 'J [Acm^{-2}]'
set xlabel 'U_{bias} [meV]'

# set yrange [0:]
# set xrange [0:10]
# set logscale y

set format y "%.1E"

set key bottom right

AU_nm = 0.0529; E0 = 1.602E-19; AU_cm2 = 2.8e-17; AU_cm = 5.29e-9

plot '../ivChar_uBarr0meV.out' u ($1*1e3):2 w l ls 1 title "U_0 = 0 meV",\
    '../ivChar_uBarr10meV.out' u ($1*1e3):2 w l ls 2 title "U_0 = 10 meV",\
    '../ivChar_uBarr20meV.out' u ($1*1e3):2 w l ls 3 title "U_0 = 20 meV"
