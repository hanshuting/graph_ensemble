#!/bin/bash
# TODO: Error check

# Submission strings for Yeti
#PBS -W group_list=yetidsi
#PBS -l mem=32000mb
#PBS -l walltime=04:00:00
#PBS -M kt2384@columbia.edu

matlab -nodesktop -nosplash -r "runDataParsing; exit;"

