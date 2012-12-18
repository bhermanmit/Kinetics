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
set title 'Constant Reactivity, Reference at 1e-5 s'
plot '-' using 1:2 with points pointtype 5 pointsize 1 linecolor rgb 'blue' title '1st Order Implicit Euler', \
     '-' using 1:2 with points pointtype 7 pointsize 1 linecolor rgb 'red'  title '4th Order Gen. Runge Kutta', \
     '-' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb 'blue' title 'Exact 1st Order', \
     '-' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb 'red' title 'Exact 4th Order'
4 7.288681e-01
2 3.580820e-01
1 1.775216e-01
5.000000e-01 8.838947e-02
2.500000e-01 4.410297e-02
1.250000e-01 2.202868e-02
e
4 3.033211e-05
2 4.355798e-06
1 3.064239e-07
5.000000e-01 2.040745e-08
2.500000e-01 1.312807e-09
1.250000e-01 7.753881e-11
e
4 7.049177e-01
2 3.524589e-01
1 1.762294e-01
5.000000e-01 8.811471e-02
2.500000e-01 4.405736e-02
1.250000e-01 2.202868e-02
e
4 8.130533e-05
2 5.081583e-06
1 3.175989e-07
5.000000e-01 1.984993e-08
2.500000e-01 1.240621e-09
1.250000e-01 7.753881e-11
e
