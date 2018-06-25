graph_ensemble
==============

by Shuting Han, Apr. 2017. Updated by Jonathan Shor, June 2018.

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
The codebase is organized as following: all experiments should be under `fwMatch-darpa/expt/` directory; `scripts/` has my matlab code to process the results; `fwMatch-darpa/src/` has all source code needed to run the model; `fwMatch-darpa/thirdparty/` has some third-party packages needed for the model.

## Working with Yeti
See documentation at https://wikis.cuit.columbia.edu/confluence/display/rcs/Yeti+HPC+Cluster+User+Documentation.

## Data format
A binary spike matrix should be stored in a `.mat` file, under the variable name `data`.
`data` needs to be number of frame by number of neuron.
The `.mat` file should be named as `<experiment>_<condition>.mat`, for example, given the `m21_d2_vis` dataset, the “experiment” is `m21_d2_vis`, the condition could be something like `high_add_neuron`, therefore the `.mat` file should be named as `m21_d2_vis_high_add_neuron.mat`.
This allows you to run multiple files that are originated from the same dataset but processed differently (high activity frame vs all frames, visual stimulations only vs all frames, no add neuron vs add neuron model, etc.) at the same time.
In the simple example below, we'll use an experiment `test` and a condition `1`.

All of the `.mat` files should be saved in the same directory, for example `~/data/<filename>` for each. For all runs under the same experiment name, this is required.

## Running CRF model - An example
1. Upload a data file `~/data/test_1.mat`
2. Go to expt directory: `cd fwMatch-darpa/expt`
3. Make a template dir for your experiment: `cp -r temporal_template/ test_template/`
This is only required once per experiment name.

4. Modify `create_script.pl` to your values:
```
$EXPT_NAME = "test";
$DATA_DIR = "~/data/";
$USER = "UNI";
$EMAIL = "UNI\@columbia.edu";
```
NOTE the backslash preceding the @ in the email address.

Between line 29 and 43, you can also adjust the model parameter test ranges here.
See below for an explanation of the parameters.

5. Start an interactive job:
```
qsub -I -q interactive -W group_list=yetibrain -l walltime=00:30:00,mem=2000mb
```
6. Run create_script, passing the condition name to it. For our example:
```
./create_script.pl 1
```
This script will create a working directory for this run based on the experiment template directory; the working directory will be named as `<experiment>_<condition>_loopy/` (in this case, `test_1_loopy/`).
Then, it will write the following files to your experiment directory: `get_real_data.m` (for loading data), `write_config_for_loopy.m` (configuration file template for the model).
The next thing it does is to start matlab, and write configuration files for each parameter combination by executing the function `create_config_files.m`.
The default settings will generate 30 files named config1.m through config30.m under the directory.
Finally, it brings you back to the `expt/` directory.

If you want to modify your YETI submission script settings (for example, letting YETI send you an email when the job is done), modify `create_config_files.m` in the template directory before running `create_script.pl`; you can safely delete the whole working directory and rerun `create_script.pl` after making edits if needed.
Between line 128 and 158 is where you look for YETI-related information.
For Luis’s datasets, 12 hours walltime and 8G memory is enough.

7. Go to working directory and start job:
```
cd test_1_loopy/
./start_jobs.sh
```
8. Check your job status by: `qstat -t -u [your_UNI]`
9. Once the job is done running, start an interactive job, and do the following:
```
matlab -nodesktop -nosplash -nodisplay
addpath(genpath(‘your/path/to/this/repo’));     % For example, '~/graph_ensemble'
cd fwMatch-darpa/expt/test_1_loopy/     % working directory
merge_all_models;
save_best_model;
```
Typing `best_model`, and then `time_span` should display the model parameters.
Take a note of `s_lambda`, `p_lambda`, `density`, and `time_span` if you want to run shuffled controls of this dataset.
A `model_collection.mat` and `<experiment>_<condition>_loopy_best_model_full.mat` file should be saved under `results/`.
Finally, exit matlab and finish the interactive job:
```
exit;
logout
```

## Running shuffled controls - An example
This produces 100 shuffled versions of the data file, and then trains a model on each using a single set of specified parameters.

1. Go to expt directory: `cd fwMatch-darpa/expt`
2. Make a template directory for your experiment name. With our example:
```
cp -r shuffled_template/ shuffled_test_template/
```
3. Modify `create_shuffle_script.pl` to your values:
```
$USER = "UNI";
$EMAIL = "UNI\@columbia.edu";
$SOURCE_DIR = "~/graph_ensemble/";
$EXPT_NAME = "test";
@EE = ("1");                        # condition name
$DATA_DIR = "~/data/";
@DENSITY = (0.29);                  # put your best density value here
@S_LAMBDA = (1.8206e-04);           # put your best s_lambda here
@P_LAMBDA = (56.2341);              # put your best p_lambda here
@TIME_SPAN = (2);                   # put your best time_span here
```
   NOTE the required backslash preceding the @ in the email address.

4. Start an interactive job, and run: `./create_shuffle_script.pl`
This will produce the working directory named as `shuffled_<experiment>_<condition>_loopy/` (in this case, `shuffled_test_1_loopy/`).
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
save_shuffled_models;
exit;
```

NOTE: If you have multiple conditions to assess, each with their own best parameters, you can process all at once by taking advantage of Perl's array notation in the following way:
```
@EE = ("1", "2", "3");    # condition names
@DENSITY = (0.29, 0.15, 0.2);  # put your best density values here
@S_LAMBDA = (1.8206e-04, 0.1, 3.5e-05); # put your best s_lambda here
@P_LAMBDA = (56.2341, 10, 1.0e2); # put your best p_lambda here
@TIME_SPAN = (2, 3, 2);
```
Notice that variables prefixed by `$` must be identical for all conditions.
Otherwise, a new script will need to be created.

## Finding core ensembles in the model
Assuming the CRF model has been trained, use `crf_core_demo.m` with the example data and model to find the ensembles correspondong to each stimulus.

## References
* [This paper - to be cited]
* Tang, K., Ruozzi, N., Belanger, D., and Jebara, T. (2016). Bethe Learning of Graphical Models via MAP Decoding. Artificial Intelligence and Statistics (AISTATS).
