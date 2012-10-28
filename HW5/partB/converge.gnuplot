set terminal pdf enhanced 
set output "partB.pdf"
set key box linestyle 1 linecolor rgb "black"
set log y
set grid
set xlabel "Iteration Number [-]"
set ylabel "L-2 Norm [-]"
plot 'converge_pj.out' using 1:2 with lines linewidth 3.0 title "Point-Jacobi", \
     'converge_gs.out' using 1:3 with lines linewidth 3.0 title "Gauss-Seidel", \
     'converge_sor.out' using 1:2 with lines linewidth 3.0 title "SOR {/Symbol w}=1.6"
