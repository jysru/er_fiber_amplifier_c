#!/usr/bin/gnuplot

reset
set terminal pdfcairo
set output "Pp.pdf"
set title "Pp"
set xlabel "z / m"
set ylabel "Pp / W"
set grid
plot "dataZ.dat" using 1:2 w l lw 3 title "Pp"


set output "N.pdf"
set title "N1 et N2"
set xlabel "z / m"
set ylabel "m-3"
plot "dataZ.dat" using 1:4 w l lw 3 title "N1", "dataZ.dat" using 1:5 w l lw 3 title "N2"


set output "sigma.pdf"
set title "N1 et N2"
set xlabel "lambda / nm"
set ylabel "m2"
plot "sigma.txt" using 1:2 w l lw 3 title "Sigma_E", "sigma.txt" using 1:3 w l lw 3 title "Sigma_A"