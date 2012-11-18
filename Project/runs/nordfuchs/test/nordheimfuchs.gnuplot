set terminal pdf dashed
set output 'nordheimfuchs.pdf'
# set terminal wxt persist
rho = 3.5*0.006648
beta = 0.006648
pnl = 1.86638609869525646E-005
K = 1.0/(148.3*312)
a = 2.5e-5
w = (rho - beta)/pnl
nhat = pnl*w**2/(2*a*K)
sech(x) = 1/cosh(x)
n(x) = nhat*(sech((w*(0.038125-x))/2))**2
T(x) = -(pnl*w)/a*tanh((w*(0.038125-x))/2) + 6.7152e+02 + 900
r(x) = beta + (rho - beta)*tanh((w*(0.038125-x))/2) - (beta + (rho - beta)*tanh((w*(0.038125))/2)) + rho
set multiplot layout 2,2
set xtics font ",5.5"
set ytics font ",5.5"
set xrange[0:0.1]
set key font ",3"
set key spacing 0.7 
set key samplen 2
set ylabel offset 2.0
set xlabel font ",5.5"
set ylabel font ",5.5"
set xlabel "Time [s]"
set ylabel "Power [W]"
plot 'data' using 1:3 with lines linewidth 3.0 title "Runge-Kutta NFM", 'data' using 1:(n($1)) with lines linewidth 3.0 linetype 2 linecolor rgb "blue" title "Anaytic NFM"
set key bottom right
set xlabel "Time [s]"
set ylabel "Fuel Temperature [K]"
plot 'data' using 1:4 with lines linewidth 3.0 title "Runge-Kutta NFM", 'data' using 1:(T($1)) with lines linewidth 3.0 linetype 2 linecolor rgb "blue" title "Analytic NFM"
set key top right
set ylabel offset 1.0
set xlabel "Time [s]"
set ylabel "Reactivity [$]"
plot 'data' using 1:2 with lines linewidth 3.0 title "Runge-Kutta NFM", 'data' using 1:(r($1)/beta) with lines linewidth 3.0 linetype 2 linecolor rgb "blue" title "Analytic NFM"
unset key
set ylabel offset 2.0
set xlabel "Time [s]"
set ylabel "Timestep [s]"
plot 'data' using 1:5 with lines linewidth 3.0
