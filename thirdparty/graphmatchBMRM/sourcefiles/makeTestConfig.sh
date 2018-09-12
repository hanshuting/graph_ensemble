#!/bin/bash
baseNames=(hotel house)
gaps=(0 10 20 30 40 50 60 70 80 90)

for baseName in ${baseNames[*]}; do
  testLosses="${baseName}_test_loss.txt"
  cat /dev/null > $testLosses

  lambdasStr=$(<"${baseName}_lambdas.txt")
  lambdas=($lambdasStr)

  for ((i=0; i<${#gaps[*]}; ++i)); do
    gap=${gaps[$i]}
    lambda=${lambdas[$i]}

    echo $baseName $gap $lambda

    configFile="${baseName}_gap${gap}_lambda${lambda}_linear.config.TEST"
    sed -e "s/{baseName}/$baseName/g" -e "s/{gap}/$gap/g" -e "s/{lambda}/$lambda/g" graphmatch_linear.config.TEST.template > $configFile

    testLog=${configFile/config/testLog}
    ./linear-bmrm-predict $configFile > $testLog

    testLearnLoss=$(grep "with weights" $testLog | cut -d ':' -f 2)
    echo $testLearnLoss >> $testLosses
  done
done

