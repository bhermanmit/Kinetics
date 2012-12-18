set terminal pdf dashed
set output 'order.pdf'
set key bottom right
set key box linestyle 1 linecolor rgb 'black'
set key width -3.0
set log x
set log y
set format y '%7.1e'
set mxtics 10
set grid x y mxtics
set xlabel 'Time Step [s]'
set ylabel 'Difference from Reference [%]'
set title 'Ramp Reactivity 0.1$ over 10s, Reference timestep: 1e-5s time: 10s'
plot '-' using 1:2 with points pointtype 7 pointsize 1 linecolor rgb 'red'  title '4th Order Gen. Runge Kutta', \
     '-' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb 'red' title 'Exact 4th Order'
2 5.920295e-02
1 1.528475e-02
5.000000e-01 3.793634e-03
2.500000e-01 9.011907e-04
1.250000e-01 2.016295e-04
1.000000e-02 3.688790e-07
1.000000e-03 1.807917e-10
1.000000e-04 3.636363e-10
e
2 2.892666e+03
1 1.807917e+02
5.000000e-01 1.129948e+01
2.500000e-01 7.062174e-01
1.250000e-01 4.413859e-02
1.000000e-02 1.807917e-06
1.000000e-03 1.807917e-10
1.000000e-04 1.807917e-14
e
