set terminal pdf dashed
set output 'pk.pdf'
set xtics font ",8"
set ytics font ",8"
set y2tics font ",8"
# set xrange[0:0.1]
set key font ",8"
# set key spacing 0.7 
# set key samplen 2
# set ylabel offset 2.0
set xlabel font ",8"
set ylabel font ",8"
set y2label font ",8"
set xlabel "Time [s]"
set ylabel "Relative Power [-]"
set ytics nomirror
set y2tics
set y2label "Reactivity [$]"
set key top left
plot 'data' using 1:3 with lines linewidth 4.0 linecolor rgb "blue" title "Power" axes x1y1, 'data' using 1:2 with lines linewidth 4.0 linecolor rgb "red" title "Reactivity" axes x1y2
