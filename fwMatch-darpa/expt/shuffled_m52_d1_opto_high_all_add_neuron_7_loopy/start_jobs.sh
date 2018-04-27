#If you are having problem with line endings use ":set ff=unix" in vim

rm -f ./results/result*.mat
rm -f ./yeti_logs/*
rm -f ./job_logs/*
cd ../.. && qsub expt/shuffled_m52_d1_opto_high_all_add_neuron_7_loopy/yeti_config.sh
