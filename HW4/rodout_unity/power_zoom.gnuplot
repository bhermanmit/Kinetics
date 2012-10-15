set terminal pdf dashed
set output "power_rodout_unity_zoom.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Rod OUT Shape Function - Unity Weighted (zoom)"
set grid
set key bottom right
set xrange[0:10]
set yrange[0.95:1.1]
plot 'kinetics_power.out' using 1 with lines linetype 1 linecolor rgb "red" linewidth 3.0 title "Reference", 'pke_power.out' using 1 with lines linetype 1 linecolor rgb "green" linewidth 3.0 title "Classic PK", 'gpke_power.out' using 1 with lines linetype 1 linecolor rgb "blue" linewidth 3.0 title "2-grp PK"
q
