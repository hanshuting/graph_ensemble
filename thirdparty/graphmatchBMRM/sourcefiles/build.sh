#!/bin/bash

[ -h "graphmatchloss.cpp" ] || ln -s ../bmrm-2.1/loss/graphmatchloss.cpp .
[ -h "graphmatchloss.hpp" ] || ln -s ../bmrm-2.1/loss/graphmatchloss.hpp .
[ -h "graphdata.cpp" ] || ln -s ../bmrm-2.1/data/graphdata.cpp .
[ -h "graphdata.hpp" ] || ln -s ../bmrm-2.1/data/graphdata.hpp .

cd ../bmrm-2.1/linear-bmrm/

make && cp linear-bmrm-predict ../../sourcefiles/ && cp linear-bmrm-train ../../sourcefiles/
