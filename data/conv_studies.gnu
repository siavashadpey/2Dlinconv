set terminal postscript portrait color enhanced font "Helvetica,12" 
set pointsize 1.5
set size square 0.8,0.8
set key bottom right
set key samplen 1

set style line 1  linetype 3 pt 4  lc rgb "#206da1" lw 1
set style line 2  linetype 4 pt 6  lc rgb "#af584b" lw 1

set style line 3  linetype 3 pt 8  lc rgb "#ff7a24" lw 1 
set style line 4  linetype 3 pt 10 lc rgb "#ffd14f" lw 1
set style line 5  linetype 3 pt 12 lc rgb "#a88ebd" lw 1
set style line 6  linetype 3 pt 14 lc rgb "#fe1832" lw 1 
set style line 7  linetype 3 pt 0  lc rgb "#010000" lw 1

set style line 8  linetype 3 pt 8  lc rgb "#3BAE5C" lw 1
set style line 9  linetype 3 pt 10 lc rgb "#3644AB" lw 1
set style line 10 linetype 3 pt 12 lc rgb "#E79D9D" lw 1
set style line 11 linetype 3 pt 14 lc rgb "#8F857D" lw 1

set style line 12 linetype 3 pt 10 lc rgb "#880069" lw 1
set style line 13 linetype 3 pt 12 lc rgb "#86FF8D" lw 1
set style line 14 linetype 3 pt 14 lc rgb "#581845" lw 1

set style line 15 linetype 3 pt 8  lc rgb "#687B5A" lw 1
set style line 16 linetype 3 pt 10 lc rgb "#00BFFF" lw 1
set style line 17 linetype 3 pt 12 lc rgb "#333300" lw 1
set style line 18 dt 2

set logscale 
set format y "10^{%L}"
set key top left
set output "conv_studies.eps"
set xlabel "h"
set ylabel "||u - u_h||_H"
plot 'conv_stud_gamma_1.dat'   using 1:2 with linespoints linestyle 1 title '{/Symbol G} p=1'   , \
     'conv_stud_gamma_2.dat'   using 1:2 with linespoints linestyle 2 title '{/Symbol G} p=2'   , \
     'conv_stud_omega_1.dat'   using 1:2 with linespoints linestyle 3 title '{/Symbol W} p=1'   , \
     'conv_stud_diagE_1.dat'   using 1:2 with linespoints linestyle 4 title '{/Symbol X} p=1'   , \
     'conv_stud_gamma_1.dat'   using 1:3 with linespoints linestyle 7 dt 2 notitle, \
     'conv_stud_gamma_2.dat'   using 1:3 with linespoints linestyle 7 dt 2 notitle, \
     'conv_stud_omega_1.dat'   using 1:3 with linespoints linestyle 7 dt 2 notitle, \
     'conv_stud_diagE_1.dat'   using 1:3 with linespoints linestyle 7 dt 2 notitle, \
     'conv_stud_gamma_1.dat'   using 1:3:4 with labels offset -2, char 1 notitle, \
     'conv_stud_gamma_2.dat'   using 1:3:4 with labels offset -2, char 1 notitle, \
     'conv_stud_omega_1.dat'   using 1:3:4 with labels offset -2, char 1 notitle, \
     'conv_stud_diagE_1.dat'   using 1:3:4 with labels offset -2, char 1 notitle