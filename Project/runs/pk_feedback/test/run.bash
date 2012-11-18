../../../src/cmfd > output
grep "POWER" output | awk {'print $2 " " $4 " " $6 " " $8 " " $10" " $12'} > data
gnuplot pkfeedback.gnuplot 
