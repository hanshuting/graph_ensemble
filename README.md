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
The original version from Kui Tang in Jebara’s group is at https://github.com/kuitang/fwmatch-public.
This repo has a modified version.
After downloading the codebase to your yeti directory, compile the dependencies as described in `fwMatch-darpa/README.md`.

The codebase is organized as following: all experiments should be under `fwMatch-darpa/expt/` directory; `scripts/` has matlab code to process the results; `fwMatch-darpa/src/` has all source code needed to run the model; `fwMatch-darpa/thirdparty/` has some third-party packages needed for the model.

## Working with Yeti
See documentation at https://wikis.cuit.columbia.edu/confluence/display/rcs/Yeti+HPC+Cluster+User+Documentation.

Ensure you are using a current version of MATLAB; version R2016b and later are known to work.
Executing
```
module load matlab/2017a
```
will ensure the R2017a version is used until you logout.

## Data format
Two binary spike matrices should be stored in a `.mat` file, under the variable names `data` and `stimuli`.
`data` needs to be number of frame by number of neuron, corresponding to the neuron recording.
`stimuli` needs to be number of frame by number of stimuli, corresponding to the when each stimuli was present.

The `.mat` file should be named as `<experiment>_<condition>.mat`, for example, given the `m21_d2_vis` dataset, the “experiment” is `m21_d2_vis`, the condition could be something like `high_add_neuron`, therefore the `.mat` file should be named as `m21_d2_vis_high_add_neuron.mat`.
This allows you to run multiple files that are originated from the same dataset but processed differently (high activity frame vs all frames, visual stimulations only vs all frames, no add neuron vs add neuron model, etc.) at the same time.
In the simple example below, we'll use an experiment `test` and a condition `1`.

All of the `.mat` files should be saved in the same directory, for example `~/data/<filename>` for each.
For all runs under the same experiment name, this is required.

## Running CRF model - An example
1. Upload a data file. Using an experiment name of "test" and a condition name of "1" as an example, we might upload to `~/data/test_1.mat`.
2. Go to expt directory. From the root directory of this repo: `cd fwMatch-darpa/expt`
3. Make a template dir for your experiment: `cp -r temporal_template/ test_template/`
This is only required once per experiment name.

4. Modify `create_script.pl` to your values:
   ```
   $EXPT_NAME = "test";
   $DATA_DIR = "~/data/";
   $USER = "UNI";
   $EMAIL = "UNI\@columbia.edu";
   ```
   NOTE the required backslash preceding the @ in the email address.

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
9. Once the job is done running, do the following from an **interactive job** (see step 5):
   ```
   matlab -nodesktop -nosplash -nodisplay
   addpath(genpath(‘your/path/to/this/repo’));     % For example, '~/graph_ensemble'
   cd fwMatch-darpa/expt/test_1_loopy/             % working directory
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

4. Start an **interactive job**, and run: `./create_shuffle_script.pl`
This will produce the working directory named as `shuffled_<experiment>_<condition>_loopy/` (in this case, `shuffled_test_1_loopy/`).
After this step, you can finish the interactive job by typing logout.

5. Go to your working directory, and run the job that generates shuffled dataset first:
   ```
   cd shuffled_test_1_loopy/
   ./shuffle_start_job.sh
   ```

   This produces 100 shuffled versions of the dataset.
   Where this occurs is defined by `create_shuffle_script.pl`, mostly the section that writes `gn_shuff_data.m` as this is the matlab script that does the actual shuffling and saving.
   The location they are saved is defined by `<data dir>/shuffled/shuffled_<experiment>_<condition>_loopy/`, so in this example it would be `~/data/shuffled/shuffled_test_1_loopy/`.  
    Settings for the yeti scheduler can be found in the section that writes `shuffle_yeti_config.sh`.

   NOTE: generated files will overwrite any existing files with identical names.

6. Monitor this job (usually it’s done within an hour).
Once it’s finished, go to working directory, and start training CRF models on the shuffled datasets:
   ```
   ./start_jobs.sh
   ```
7. When all jobs are done, start an **interactive job**, and start matlab:
   ```
   matlab -nodesktop -nosplash -nodisplay
   addpath(genpath(‘your/path/to/this/repo’));              % For example, '~/graph_ensemble'
   cd fwMatch-darpa/expt/shuffled_test_1_loopy/             % working directory
   merge_all_models;
   save_shuffled_models;
   exit;
   ```

   This will compile the trained models into two files: `model_collection.mat` and `shuffled_<experiment>_<condition>_loopy_fulldata.mat`, both in the `results/` folder in the working directory.

NOTE: You can simultaneously process multiple conditions, each with their own best parameters, by taking advantage of Perl's array notation when editing `create_shuffle_script.pl` in step 3.
Each condition will have its own working directory, config files, etc., as normal.
Here is an example for three conditions named "1", "2", and "3":
```
@EE = ("1", "2", "3");                  # condition names
@DENSITY = (0.29, 0.15, 0.2);           # put your best density values here
@S_LAMBDA = (1.8206e-04, 0.1, 3.5e-05); # put your best s_lambdas here
@P_LAMBDA = (56.2341, 10, 1.0e2);       # put your best p_lambdas here
@TIME_SPAN = (2, 3, 2);                 # put your best time_spans here
```

Parameters should be ordered by condition consistently for each array variable assignment.  
Variables prefixed by `$` must be identical for all conditions, unlike `@` prefixed array variables.
Otherwise, a new script will need to be created.

## Finding core ensembles in the model
With a trained CRF model on the dataset of interest, and a collection of models trained on shuffled versions of the dataset, use `scripts/core/find_temporal_ens_nodes.m` to find the ensembles corresponding to each stimulus.
We will continue our previous example:

1. Start an interactive job:
   ```
   qsub -I -q interactive -W group_list=yetibrain -l walltime=00:30:00,mem=2000mb
   ```
2. Start matlab and load the models and data:
   ```
   matlab -nodesktop -nosplash -nodisplay
   addpath(genpath(‘your/path/to/this/repo’));             % For example, '~/graph_ensemble'
   cd fwMatch-darpa/expt/
   best_model = load('<experiment>_<condition>_loopy/results/<experiment>_<condition>_loopy_best_model_full.mat');
   shuffle_model = load('shuffled_<experiment>_<condition>_loopy/results/shuffled_<experiment>_<condition>_loopy_fulldata.mat');
   load('~/data/<experiment>_<condition>.mat');            % Loads variables `data` and `stimuli`
   ```
3. Find ensemble nodes:
   ```
   ens_nodes = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli)
   ```
   `ens_nodes` is a cell vector where each cell contains the ensemble neurons found for each stimuli.
   Each such stimuli cell contains a further cell vector where each cell contains the ensemble neurons found for each offset frame of the `time_span` window.

Another script, `scripts/core/find_temporal_crf_core.m`, can also be used to find ensemble neurons and plot some features, including spatial arrangement if coordinates are provided.

## References
* [This paper - to be cited]
* Tang, K., Ruozzi, N., Belanger, D., and Jebara, T. (2016). Bethe Learning of Graphical Models via MAP Decoding. Artificial Intelligence and Statistics (AISTATS).
