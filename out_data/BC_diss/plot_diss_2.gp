set terminal pngcairo dashed enhanced size 800, 500
set output 'diss_BC.png'

set multiplot

set style line 1 lc rgb '#000000' lw 2 dashtype 1
set style line 2 lc rgb '#000080' lw 2 dashtype 1
set style line 3 lc rgb '#FF0000' lw 2 dashtype 1
set style line 4 lc rgb '#008000' lw 2 dashtype 1
set style line 5 lc rgb '#800000' lw 2 dashtype 1
set style line 6 lc rgb '#800080' lw 2 dashtype 1
set style line 7 lc rgb '#008080' lw 2 dashtype 1
set style line 8 lc rgb '#808000' lw 2 dashtype 1

#
# Plot 1
#

set size .99, 1       # set size of inset
set origin 0, 0

set ylabel 'J({/Symbol t})/J({/Symbol t}=0)'  # 'I({/Symbol t})/J_0'
set xlabel '{/Symbol t}_{/Symbol G} [s]'

# set yrange [0:1.8]
# set xrange [1e-15:1e-10]

# set key title "L[f; {/Symbol t}]"  # '{/Symbol t}'
set key title 'Profil'
set key top left

# set title 'Transport characteristics'

set logscale x

AU_nm = 0.0529; E0 = 1.602E-19; AU_Acm2 = 2.364e14; AU_cm2 = 2.8e-17
file = ''

plot 'diss_lorentz.out' u 2:3 w l ls 2 title 'Lorentz',\
    'diss_gauss.out' u 2:3 w l ls 3 title 'Gauss',\
    'diss_voigt.out' u 2:3 w l ls 4 title 'Voigt'

#
# Inset
#

set size 0.56, 0.56       # set size of inset
set origin 0.35, 0.12

set ylabel 'N({/Symbol t})/N({/Symbol t}=0)'
# set yrange [0.8:2]
# set xrange [1e-14:1e-10]

unset key

plot 'diss_lorentz.out' u 2:4 w l ls 2 title 'Lorentz',\
    'diss_gauss.out' u 2:4 w l ls 3 title 'Gauss',\
    'diss_voigt.out' u 2:4 w l ls 4 title 'Voigt'
