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
`data` needs to be number of frames by number of neurons, corresponding to the neuron recordings.
`stimuli` needs to be number of frames by number of stimuli, corresponding to when each stimuli is present.
Submitting a stimulus with no occurrences is not accepted, and needs to be reoved from the matrix.

The `.mat` file should be named as `<experiment>_<condition>.mat`, for example, given the `m21_d2_vis` dataset, the “experiment” is `m21_d2_vis`, the condition could be something like `high_add_neuron`, therefore the `.mat` file should be named as `m21_d2_vis_high_add_neuron.mat`.
This allows you to run multiple files that are originated from the same dataset but processed differently (high activity frame vs all frames, visual stimulations only vs all frames, no add neuron vs add neuron model, etc.) at the same time.
In the simple example below, we'll use an experiment `test` and a condition `1`.

All of the `.mat` files should be saved in the same directory, for example `~/data/<filename>` for each.
For all runs under the same experiment name, this is required.

## Running CRF model - An example
1. Upload a data file. Using an experiment name of "test" and a condition name of "1" as an example, we might upload to `~/data/test_1.mat`.
2. Go to expt directory. From the root directory of this repo: `cd fwMatch-darpa/expt`
3. Edit the "USER EDITABLE VARIABLES" in `run_full_temporal.py` to your values.
   `SOURCE_DIR` should refer to the root directory of where this repo is installed; typically this directory is named graph_ensemble.
   For this example, we might use the following values instead of the defaults:
   ```
   EXPT_NAME = "test"
   USER = "UNI"
   EMAIL = "UNI@columbia.edu"
   ```
4. Run `run_full_temporal` with Python 3.5 or greater, passing the condition name to it. For our example with a condition name of 1:
   ```
   python3 run_full_temporal.py 1
   ```

This script will conduct a grid search across the parameter ranges specified, training a CRF model on the data file for each parameter combination.
A working directory will be created and named as `<experiment>_<condition>_loopy/` (in this case, `test_1_loopy/`), and a directory within it named `results` will contain the trained models and the best model.

The best parameters will be extracted and used to produce shuffled control datasets, the number of which is controlled by the NSHUFFLE in `run_full_temporal.py`.
Another working directory for the shuffle models will be created and named as `shuffled_<experiment>_<condition>_loopy/` (in this example, `shuffled_test_1_loopy/`).
Again, a `results` directory inside will contain the trained models.


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
   best_model = load('<experiment>_<condition>_loopy/results/best_model_full.mat');
   shuffle_model = load('shuffled_<experiment>_<condition>_loopy/results/fulldata.mat');
   load('~/data/<experiment>_<condition>.mat');            % Loads variables `data` and `stimuli`
   ```
3. Find ensemble nodes:
   ```
   ens_nodes = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli)
   ```
   `ens_nodes` is a cell vector where each cell contains the ensemble neurons found for each stimuli.
   Each such stimuli cell contains a further cell vector where each cell contains the ensemble neurons found for each offset frame of the `time_span` window.

Another script, `scripts/core/find_plot_temporal_crf_core.m`, can also be used to find ensemble neurons and plot some features, including spatial arrangement if coordinates are provided.

## References
* Carrillo-Reid, L.\*, Han, S.\*, Taralova, E., Jebara, T., Yuste, R. (2017). Identification and Targeting of Cortical Ensembles. bioRxiv. doi: https://doi.org/10.1101/226514
* Tang, K., Ruozzi, N., Belanger, D., and Jebara, T. (2016). Bethe Learning of Graphical Models via MAP Decoding. Artificial Intelligence and Statistics (AISTATS).
