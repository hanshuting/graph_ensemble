#!/bin/bash

for file in *config; do
  trainLog=${file/config/trainLog}
  validLog=${file/config/validLog}
  lossLog=${file/config/lossLog}
  echo $trainLog
  echo $validLog

  ./linear-bmrm-train $file   > $trainLog
  ./linear-bmrm-predict $file > $validLog

  nIters=$(grep iterations $trainLog | cut -d ':' -f 2)
  learnTime=$(grep Total $trainLog | cut -d ':' -f 2)
  trainLearnLoss=$(grep "with weights" $trainLog | cut -d ':' -f 2)

  validLearnLoss=$(grep "with weights" $validLog | cut -d ':' -f 2)

  echo $nIters $learnTime $trainLearnLoss $validLearnLoss > $lossLog

done

