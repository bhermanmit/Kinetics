set terminal pdf dashed
set output "error_exact_unity.pdf"
set key bottom right
set xlabel "Time [s]"
set ylabel "Relative Core Power [-]"
set title "Comparison of Classic PKEs with PNL(t)"
set grid
plot '<paste kinetics_power.out old_l/pke_oldl.out' using ($0/50):(($2-$1)/$1) with lines linetype 1 linecolor rgb "red" linewidth 3.0 title "Old Definition - PNL(-)", \
     '<paste kinetics_power.out old_lt/pke_oldlt.out' using ($0/50):(($2-$1)/$1) with lines linetype 1 linecolor rgb "green" linewidth 3.0 title "Old Definition - PNL(t)", \
     '<paste kinetics_power.out new_l/pke_newl.out' using ($0/50):(($2-$1)/$1) with lines linetype 2 linecolor rgb "blue" linewidth 3.0 title "New Definition - PNL(-)", \
     '<paste kinetics_power.out new_lt/pke_newlt.out' using ($0/50):(($2-$1)/$1) with lines linetype 3 linecolor rgb "black" linewidth 3.0 title "New Definition - PNL(t)"
q
