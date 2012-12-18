set terminal pdf dashed
set output 'order.pdf'
set key bottom right
set key box linestyle 1 linecolor rgb 'black'
set log x
set log y
set format y '%7.1e'
set mxtics 10
set grid x y mxtics
set xlabel 'Time Step [s]'
set ylabel 'Difference from Reference [%]'
set title 'PK with feedback, Reference at time step: 1e-5s'
plot '-' using 1:2 with points pointtype 7 pointsize 1 linecolor rgb 'red'  title 'GRK4 @ t=0.002s', \
     '-' using 1:2 with lines linetype 1 linewidth 2 linecolor rgb 'red' title 'Exact 4th Order', \
     '-' using 1:2 with points pointtype 9 pointsize 1 linecolor rgb 'blue' title 'GRK4 @ t=2s', \
     '-' using 1:2 with points pointtype 11 pointsize 1 linecolor rgb 'orange' title 'GRK4 @ t=0.09s'
2.000000e-03 1.954779e-03
1.000000e-03 1.689959e-04
5.000000e-04 1.276439e-05
2.500000e-04 8.834061e-07
1.250000e-04 5.781205e-08
e
2.000000e-03 3.788770e-03
1.000000e-03 2.367981e-04
5.000000e-04 1.479988e-05
2.500000e-04 9.249927e-07
1.250000e-04 5.781205e-08
e
2.000000e-03 1.7813e-03
1.000000e-03 1.3043e-03
5.000000e-04 8.7778e-04
2.500000e-04 4.9802e-04
1.250000e-04 2.5607e-04
e
2.000000e-03 3.9512e-02
1.000000e-03 2.1651e-02
5.000000e-04 1.1291e-02
2.500000e-04 5.6840e-03
1.250000e-04 2.7621e-03
e

