set terminal pdf dashed
set output 'pkfeedback.pdf'
# set terminal wxt persist
set multiplot layout 2,2
set xtics font ",5.5"
set ytics font ",5.5"
set xrange[0:2]
unset key
set format y "%7.1e"
set ylabel offset 2.0
set xlabel font ",5.5"
set ylabel font ",5.5"
set xlabel "Time [s]"
set ylabel "Power [W]"
plot 'data' using 1:3 with lines linewidth 3.0
set key center right
set key font ",4"
set key samplen 2
set xlabel "Time [s]"
set ylabel "Temperature [K]"
plot 'data' using 1:4 with lines linewidth 3.0 title "Fuel", 'data' using 1:5 with lines linewidth 3.0 linecolor rgb 'blue' title "Coolant"
unset key
set ylabel offset 1.0
set xlabel "Time [s]"
set ylabel "Reactivity [$]"
plot 'data' using 1:2 with lines linewidth 3.0
set ylabel offset 2.0
set xlabel "Time [s]"
set ylabel "Timestep [s]"
plot 'data' using 1:6 with lines linewidth 3.0
