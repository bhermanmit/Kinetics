set terminal pdf
set output "power_exact_unity.pdf"
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Comparison of Classic PKEs with PNL(t)"
set grid
plot 'kinetics_power.out' using ($0/50):1 with lines linewidth 3.0 title "Reference", \
     'old_l/pke_oldl.out' using ($0/50):1 with lines linewidth 3.0 title "Old Definition - PNL(-)", \
     'old_lt/pke_oldlt.out' using ($0/50):1 with lines linewidth 3.0 title "Old Definition - PNL(t)", \
     'new_l/pke_newl.out' using ($0/50):1 with lines linewidth 3.0 title "New Definition - PNL(-)", \
     'new_lt/pke_newlt.out' using ($0/50):1 with lines linewidth 3.0 title "New Definition - PNL(t)"
q
