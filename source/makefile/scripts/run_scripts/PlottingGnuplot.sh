gnuplot
plot "Mass.txt" using 1:2 title 'WM' with linespoints,"Mass.txt" using 1:3 title 'GM' with linespoints, "Mass.txt" using 1:4 title 'CSF' with linespoints
set grid
set xlabel "time" ; set ylabel "mass"; set title "Mass Conservation"
