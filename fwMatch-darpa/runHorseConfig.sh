num=$1

/exp/comm/matlab/bin/matlab -nodisplay  -r "addpath(genpath('.')); runExpt('horseRhoCrossVal',$num);"
