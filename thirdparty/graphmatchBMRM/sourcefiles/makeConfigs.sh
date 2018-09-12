#!/bin/bash
baseNames=(hotel house)
gaps=(0 10 20 30 40 50 60 70 80 90)
lambdas=(0.0001 0.0003 0.001 0.003 0.01 0.03 0.1 0.3 1 3 10 30 100)

for baseName in ${baseNames[*]}; do
  for gap in ${gaps[*]}; do
    for lambda in ${lambdas[*]}; do
#      trainFile="../Data/${baseName}/pairs/${baseName}s${gap}_train.txt"

      echo $baseName $gap $lambda
      outFile="${baseName}_gap${gap}_lambda${lambda}_linear.config"
      sed -e "s/{baseName}/$baseName/g" -e "s/{gap}/$gap/g" -e "s/{lambda}/$lambda/g" graphmatch_linear.config.template > $outFile
#      sed -e "s/{baseName}/$baseName/g" -e "s/{gap}/$gap/g" graphmatch_linear.config.template > $outFile
    done
  done
done

