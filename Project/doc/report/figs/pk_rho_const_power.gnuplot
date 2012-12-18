set terminal pdf dashed
set output 'pk.pdf'
set xlabel "Time [s]"
set ylabel "Relative Power [-]"
set ytics nomirror
set y2tics
set y2label "Reactivity [$]"
set key top left 
plot 'data' using 1:3 with lines linewidth 4.0 linecolor rgb "blue" title "Power" axes x1y1, 'data' using 1:2 with lines linewidth 4.0 linecolor rgb "red" title "Reactivity" axes x1y2
