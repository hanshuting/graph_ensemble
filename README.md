graph_ensemble
==============

by Shuting Han, Apr. 2017

Licensing Info
--------------
[To be added]
Based on code written by Tang et al. (2016) available at https://github.com/kuitang/fwmatch-public.

Overview
--------
This project applies conditional random field models in _in vivo_ two-photon calcium imaging data in mouse V1 cortex. We built CRFs on such imaging data to identify representative neuronal ensembles corresponding to specific external stimulus, and to identify core neurons that are able to perform pattern completion.
This code is intended for parallel processings on linux servers. The following instructions are designed for Columbia Yeti HPC usage.

## CRF codebase
The original version from Kui Tang in Jebara’s group is at https://github.com/kuitang/fwmatch-public. This repo has a modified version. After downloading the codebase to your yeti directory, compile the dependencies as described in the readme file.
The codebase is organized as following: all experiments should be under `expt/ directory`; `mvc/` has my matlab code to process the results; `src/` has all source code needed to run the model; `thirdparty/` has some third-party packages needed for the model.

## Working with Yeti
See documentation at https://wikis.cuit.columbia.edu/confluence/display/rcs/Yeti+HPC+Cluster+User+Documentation. 

## Data format
A binary spike matrix should be stored in a `.mat` file, under the variable name `data`. `data` needs to be number of frame by number of neuron. The `.mat` file should be named as `experiment_condition.mat`, for example, given the `m21_d2_vis` dataset, the “experiment” is `m21_d2_vis`, the condition could be something like `high_add_neuron`, therefore the `.mat` file should be named as `m21_d2_vis_high_add_neuron.mat`. This allows you to run multiple files that are originated from the same dataset but processed differently (high activity frame vs all frames, visual stimulations only vs all frames, no add neuron vs add neuron model, etc.) at the same time. In my case, the data file is saved in `~/data/luis/[filename]`.

## Running CRF model - An example
1. Upload a data file `$DATAPATH/test_1.mat` [To be completed]
2. Go to expt directory: `cd ~/src/fwMatch-darpa/expt`
3. Make a template dir: `cp -r m21_d2_vis_template/ test_template/`
Most files in this template directory don’t matter; they’re here only for historical reasons. The only files that matter are `create_config_files.m` and `merge_all_models.m`. I will clean up the codebase in the future.
4. Modify `create_script.pl` accordingly:
```
$EXPT_NAME = "test";
@EE = ("1"); # note you can put multiple experiments here as (“1”,”2”,”3”);`
$MODEL_TYPE = "loopy";
$DATA_DIR = "/vega/brain/users/sh3276/data/luis";
```
Between line 42 and 59, modify the email address to yours, and yeti_user to your UNI. You can also adjust the model parameter test ranges here.
5. Start an interactive job:
```
qsub -I -q interactive -W group_list=yetibrain -l walltime=00:30:00,mem=2000mb
```
6. Run create_script by: ./create_script.pl
This script will create a working directory based on the template directory; the working directory will be named as `experiment_condition_loopy/` (in this case, `test_1_loopy/`). Then, it will write the following files to your experiment directory: `get_real_data.m` (for loading data), `write_config_for_loopy.m` (configuration file template for loopy models), `write_config_for_tree.m` (configuration file template for tree models). The next thing it does is to start matlab, and write configuration files for each parameter combination by executing the function `create_config_files.m`. This will generate 30 files like config1.m under the directory. Finally, it brings you back to the expt/ dir.
If you want to modify your YETI submission script settings (for example, letting YETI sending you an email when the job is done), modify `create_config_files.m` in the template directory before running `create_script.pl`. Between line 115 and 143 is where you look for YETI-related information. For Luis’s datasets, 12 hours walltime and 8G memory is enough.
7. Go to working dir and start job:
```
cd test_1_loopy/
./start_jobs.sh
```
8. Check your job status by: qstat -t -u [your_UNI]
9. Once the job is done running, start an interactive job, and do the following:
```
matlab -nodesktop -nosplash -nodisplay
addpath(genpath(‘~/src/fwMatch-darpa’));
cd ~/expt/test_1_loopy/
merge_all_models;
save_best_model;
```
Typing `best_model.theta` should display the model parameters. Take a note of `s_lambda`, `p_lambda` and density if you want to run shuffled controls of this dataset.
A model_collection.mat file should be saved under results/.
Finally, exit matlab and finish the interactive job:
```
exit;
logout
```
Before running `save_best_model` for the first time, please change the save path to your desired path. This file is located at `fwMatch-darpa/mvc/scripts/save_best_model.m`. Then a file named `test_1_loopy_best_model_full.mat` should be saved this path. 

## Running shuffled controls - An example
1. Go to expt directory: `cd ~/src/fwMatch-darpa/expt`
2. Make a template dir: `cp -r shuffled_m21_d2_vis_template/ shuffled_test_template/`
3. Modify `create_shuffle_script.pl` accordingly:
```
$EXPT_NAME = "test";
@EE = ("1");
$MODEL_TYPE = "loopy";
$DATA_DIR = "/vega/brain/users/sh3276/data/luis";
@DENSITY = (0.29);  # put your best density value here
@S_LAMBDA = (1.8206e-04); # put your best s_lambda here
@P_LAMBDA = (56.2341); # put your best p_lambda here
$NSHUFFLE = 100;
```
4. Start an interactive job, and run: `./create_shuffle_script.pl`
After this step, you can finish the interactive job by typing logout.
5. Go to your working directory, and run the job that generates shuffled dataset first:
```
cd shuffled_test_1_loopy/
./shuffle_start_job.sh
```
6. Monitor this job (usually it’s done within an hour), once it’s finished, go to working directory, and start running CRF models on shuffled data:
```
./start_job.sh
```
7. When all jobs are done, start an interactive job, and start matlab:
```
matlab -nodesktop -nosplash -nodisplay
addpath(genpath(‘~/src/fwMatch-darpa’));
cd ~/expt/shuffled_test_1_loopy/
merge_all_models;
save_shuffled_model;
exit;
```
Again, please change the save path in `save_shuffled_model.m` to your desired path before running it for the first time.

## Finding core ensembles in the model
Assuming the CRF model has been trained, use `crf_core_demo.m` with the example data and model to find the ensembles correspondong to each stimulus.

## References
* [This paper - to be cited]
* Tang, K., Ruozzi, N., Belanger, D., and Jebara, T. (2016). Bethe Learning of Graphical Models via MAP Decoding. Artificial Intelligence and Statistics (AISTATS).
