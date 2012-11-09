grep "POWER" $1 | awk '{print $2}' > "power""$2"
grep "POWER" $1 | awk '{print $4}' > "time""$2"
