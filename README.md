graph_ensemble
==============

by Shuting Han, Apr. 2017, Sep. 2018; Jonathan Shor, June 2018.

Overview
--------
This project applies conditional random field models in _in vivo_ two-photon calcium imaging data in mouse V1 cortex. We built CRFs on such imaging data to identify representative neuronal ensembles corresponding to specific external stimulus, and to identify core neurons that are able to perform pattern completion. This codebase is based on the [fwMatch](https://github.com/kuitang/fwmatch-public) repo from Jebara group.

This code was tested only on linux/macOS platforms.  
MATLAB version R2016b or later is required.  
Python 3.5+ is required for the example workflow scripts.

## Compile dependencies
After cloning the repo, run `make` in the base directory. 

Note: Please change the MEX path in `thirdparty/QPBO-v1.32.src/Makefile` (line 19) to the Matlab path in your system.
On MacOS, XCode is required for compilation.

## Data format
Two binary spike matrices should be stored in a `.mat` file, under the variable names `data` and `stimuli`.
`data` needs to be number of frames by number of neurons, corresponding to the neuron recordings.
`stimuli` needs to be number of frames by number of stimuli, corresponding to when each stimuli is present.
Submitting a stimulus with no occurrences is not accepted, and needs to be reoved from the matrix.

The `.mat` file should be named as `<experiment>_<condition>.mat`, for example, given the `m21_d2_vis` dataset, the “experiment” is `m21_d2_vis`, the condition could be something like `high_add_neuron`, therefore the `.mat` file should be named as `m21_d2_vis_high_add_neuron.mat`.
This allows you to run multiple files that are originated from the same dataset but processed differently (high activity frame vs all frames, visual stimulations only vs all frames, no add neuron vs add neuron model, etc.) at the same time, and each will have their own working directory.

## Running CRF model - An example
1. Upload a data file.
   Using an experiment name of `experiment` and a condition name of `demo` as an example (file provided with the repo), we might upload to `~/data/experiment_demo.mat`.
2. Go to expt directory. From the root directory of this repo: `cd expt`
3. Edit the "GeneralOptions" in `crf_parameters.ini` to your values.
   Ensure `experiment_name`, `data_directory`, and `source_directory` are correctly updated, and set `time_span` as desired.
   The other parameter defaults are likely reasonable for initial exploratory runs.
   For this example, we might have the following values:
   ```
   experiment_name = experiment
   data_directory = ~/data/
   source_directory = ~/graph_ensemble/
   ```
   Note that `crf_parameters.ini` is the default settings filename, but the included workflow scripts accept custom specified files.
   When using from the command line, as explained in these examples, simply pass the settings filename as the second argument after the condition name.

   **If running on Columbia's yeti cluster, see further instructions in dedicated section below.**

4. Run `run_full_temporal` with Python 3.5 or greater, passing the condition name to it. For our example with a condition name of demo:
   ```
   python3 run_full_temporal.py demo crf_parameters.ini
   ```

This script will conduct a grid search across the parameter ranges specified, training a CRF model on the data file for each parameter combination.
A working directory will be created and named as `<experiment>_<condition>/` (in this case, `experiment_demo/`), and a directory within it named `results` will contain the trained models and the best model.

The best parameters will be extracted and used to produce shuffled control datasets, the number of which is controlled by the num_shuffle option in `crf_parameters.ini`.
Another working directory for the shuffle models will be created and named as `shuffled_<experiment>_<condition>/` (in this example, `shuffled_experiment_demo/`).
Again, a `results` directory inside will contain the trained models.


## Finding core ensembles in the model
With a trained CRF model on the dataset of interest, and a collection of models trained on shuffled versions of the dataset, use `src/core/find_temporal_ens_nodes.m` to find the ensembles corresponding to each stimulus.
We will continue our previous example:

1. Start matlab. You can do so from the terminal with the following command:
   ```
   matlab -nodesktop -nosplash -nodisplay
   ```
2. Load the models and data in matlab:
   ```
   addpath(genpath(‘~/graph_ensemble/’));       % This should be your source_directory
   cd expt/
   best_model = load('experiment_demo/results/best_model_full.mat');
   shuffle_model = load('shuffled_experiment_demo/results/fulldata.mat');
   load('~/data/experiment_demo.mat');            % Loads variables `data` and `stimuli`
   ```
3. Find ensemble nodes:
   ```
   ens_nodes = find_temporal_ens_nodes(best_model, shuffle_model, data, stimuli)
   ```
   `ens_nodes` is a cell vector where each cell contains the ensemble neurons found for each stimuli.
   Each such stimuli cell contains a further cell vector where each cell contains the ensemble neurons found for each offset frame of the `time_span` window, with the first corresponding to no offset.

Another script, `scripts/core/find_plot_temporal_crf_core.m`, can also be used on a desktop system to find ensemble neurons and plot some features, including spatial arrangement if coordinates are provided.


## Working with Yeti (for Columbia users)
See yeti cluster documentation at https://wikis.cuit.columbia.edu/confluence/display/rcs/Yeti+HPC+Cluster+User+Documentation.

Ensure you are using acceptable versions of MATLAB (at least 2016b) and Python (at least 3.5).
We strongly recommend adding something like the following two lines to the `.bash_profile` file in your home directory:
```
module load matlab/2017a
module load anaconda/4.1.1-python-3.5.2
```
This will ensure the correct versions are always used.

There are several yeti specific settings sections in `crf_parameters.ini`, each named with a "Yeti" prefix, to control job submission, resource requesting, and email notifications.  
Be sure `username`, `group_id`, and `email` are updated to valid values.
The remaining defaults are likely reasonable for initial exploratory runs.

Additionally, there are yeti versions of all executable workflow scripts in the expt folder, identifiable by a `yeti_` prefix.
For example, `run_full_temporal.py` is updated to take advantage of the yeti cluster's features in `yeti_run_full_temporal.py`.
Beyond attending to the yeti settings in the settings file, execution is identical.

Note that, in the case of using ssh to connect to the cluster, being disconnected for any reason will terminate any local jobs.
Running one of the python workflow scripts directly from the command line as described above is an example of such a local job.
Consider either submitting a job to the cluster to run the script, or use a terminal multiplexer like [screen](https://linuxize.com/post/how-to-use-linux-screen/) which is already installed on the yeti cluster.


## References
* Carrillo-Reid, L.\*, Han, S.\*, Taralova, E., Jebara, T., Yuste, R. (2017). Identification and Targeting of Cortical Ensembles. bioRxiv. doi: https://doi.org/10.1101/226514
* Tang, K., Ruozzi, N., Belanger, D., and Jebara, T. (2016). Bethe Learning of Graphical Models via MAP Decoding. Artificial Intelligence and Statistics (AISTATS).
