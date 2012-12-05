grep "POWER" $1 | awk '{print $6}' > "power""$2"
grep "POWER" $1 | awk '{print $2}' > "time""$2"
