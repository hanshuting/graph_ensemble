#If you are having problem with line endings use ":set ff=unix" in vim

rm -f ./results/result*.mat
rm -f ./yeti_logs/*
rm -f ./job_logs/*
cd ../.. && qsub expt/m21_d2_vis_all_high_add_neuron_loopy/yeti_config.sh
