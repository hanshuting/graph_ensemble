for log in $(grep -l 'program exceeded' *lossLog); do
  cut -d '!' -f 2 < $log > $log.tmp
  mv $log.tmp $log
done

