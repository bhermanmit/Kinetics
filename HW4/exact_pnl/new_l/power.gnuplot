set terminal pdf
set output "power_exact_unity.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Exact PKE Parameters - Unity Weighted"
set grid
plot 'kinetics_power.out' using ($0/50):1 with lines linewidth 3.0 title "Reference", 'pke_power.out' using ($0/50):1 with lines linewidth 3.0 title "Classic PK"
q
