#!/bin/sh
#yeti_config.sh

#Torque script to run Matlab program

#Torque directives
#PBS -N m21_d2_vis_all_add_neuron_loopy
#PBS -W group_list=yetibrain
#PBS -l nodes=1:ppn=2,walltime=12:00:00,mem=8000mb
#PBS -V
#PBS -t 1-30

#set output and error directories (SSCC example here)
#PBS -o localhost:/vega/brain/users/sh3276/src/fwMatch-darpa/expt/m21_d2_vis_all_add_neuron_loopy/yeti_logs/
#PBS -e localhost:/vega/brain/users/sh3276/src/fwMatch-darpa/expt/m21_d2_vis_all_add_neuron_loopy/yeti_logs/

#Command below is to execute Matlab code for Job Array (Example 4) so that each part writes own output
cd /vega/brain/users/sh3276/src/fwMatch-darpa/
./run.sh m21_d2_vis_all_add_neuron_loopy $PBS_ARRAYID > expt/m21_d2_vis_all_add_neuron_loopy/job_logs/matoutfile.$PBS_ARRAYID
#End of script
