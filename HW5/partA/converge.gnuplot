set terminal pdf 
set output "partA.pdf"
set key box linestyle 1 linecolor rgb "black"
set key width -4.0
set log y
set grid
set xlabel "Iteration Number [-]"
set ylabel "L-2 Norm [-]"
plot 'converge.out' using 1:2 with lines linewidth 3.0 title "True Fractional Error", \
     'converge.out' using 1:3 with lines linewidth 3.0 title "Normalized Residual Error", \
     'converge.out' using 1:4 with lines linewidth 3.0 title "Successive Fractional Error"
