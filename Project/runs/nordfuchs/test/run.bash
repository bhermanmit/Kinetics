../../../src/cmfd > output
grep "POWER" output | awk {'print $2 " " $4 " " $6 " " $8 " " $10'} > data
gnuplot nordheimfuchs.gnuplot 
