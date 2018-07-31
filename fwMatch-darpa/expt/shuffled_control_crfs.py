#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import stat
import subprocess

import crf_util
from gridsearch_train_crfs import get_best_parameters

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))

# *** START USER EDITABLE VARIABLES ***
logger.setLevel(logging.DEBUG)
EXPT_NAME = "temporal"
DATA_DIR = "~/data/"
SOURCE_DIR = "~/graph_ensemble/"    # TODO: Autodetect this
USER = "jds2270"
EMAIL = "jds2270@columbia.edu"
NSHUFFLE = 100
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"
EXPT_DIR = os.path.join(SOURCE_DIR, "fwMatch-darpa", "expt")

# These parameters and their order must match best_parameters.txt.
# See save_best_parameters.m for best_parameters.txt creation.
PARAMS_TO_EXTRACT = ['s_lambda', 'density', 'p_lambda', 'time_span']
# *** end constants ***


def get_conditions_metadata(conditions):
    for name in conditions:
        experiment = "{}_{}_{}".format(EXPT_NAME, name, MODEL_TYPE)
        metadata = {'data_file': "{}_{}.mat".format(EXPT_NAME, name),
                    'experiment': experiment,
                    'shuffle_save_name': "shuffled_{}_{}".format(EXPT_NAME, name),
                    'shuffle_save_dir': os.path.join(DATA_DIR, "shuffled", experiment),
                    'shuffle_experiment': "shuffled_{}".format(experiment)
                    }
        conditions[name].update(metadata)
    return conditions


def write_shuffled_data_generating_script(experiment, data_file, save_dir, save_name):
    # script for generate shuffled data
    filepath = os.path.join(experiment, "gn_shuff_data.m")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("if exist('{}{}')~=7\n".format(DATA_DIR, save_name))
        f.write("    mkdir('{}{}');\n".format(DATA_DIR, save_name))
        f.write("end\n")
        f.write("addpath(genpath('{}'));\n".format(SOURCE_DIR))
        f.write("load(['{}{}']);\n".format(DATA_DIR, data_file))
        f.write("fprintf('Loaded: %s\\n', ['{}{}']);\n".format(DATA_DIR, data_file))
        f.write("if exist('stimuli', 'var') ~= 1\n")
        f.write("    stimuli = [];\n")
        f.write("end\n")
        f.write("num_stimuli = size(stimuli, 2)\n")
        f.write("data_raw = [data stimuli]';\n")
        f.write("for i = 1:{}\n".format(NSHUFFLE))
        f.write("\tdata = shuffle(data_raw,'exchange')';\n")
        f.write("\tstimuli = data(:, end - num_stimuli + 1:end);\n")
        f.write("\tdata = data(:, 1:end - num_stimuli);\n")
        f.write("\tsave(['{}_' num2str(i) '.mat'],'data','stimuli');\n".format(
            os.path.join(save_dir, save_name)))
        f.write("end\n")
        f.write("fprintf('done shuffling data\\n');\n")
    f.closed
    logger.info("done writing {}".format(filepath))


def write_shuffling_yeti_script(experiment):
    # write yeti script
    filepath = os.path.join(experiment, "shuffle_yeti_config.sh")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#shuffle_yeti_config.sh\n")
        f.write("#PBS -N {}\n".format(experiment))
        f.write("#PBS -W group_list=yetibrain\n")
        f.write("#PBS -l nodes=1:ppn=1,walltime=02:00:00,mem=4000mb\n")
        f.write("#PBS -m ae\n")
        f.write("#PBS -V\n")
        f.write("#set output and error directories (SSCC example here)\n")
        log_folder = "{}fwMatch-darpa/expt/{}/yeti_logs/".format(SOURCE_DIR, experiment)
        os.makedirs(os.path.expanduser(log_folder), exist_ok=True)
        f.write("#PBS -o localhost:{}\n".format(log_folder))
        f.write("#PBS -e localhost:{}\n".format(log_folder))
        f.write("#Command below is to execute Matlab code for Job Array (Example 4) " +
                "so that each part writes own output\n")
        f.write("matlab -nodesktop -nodisplay -r \"dbclear all;" +
                " addpath('{}');".format(os.path.join(EXPT_DIR, experiment)) +
                "gn_shuff_data; exit\"\n")
        f.write("#End of script\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}".format(filepath))


def write_shuffling_submit_script(experiment):
    # write submit script
    filepath = os.path.join(experiment, "shuffle_start_job.sh")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("qsub shuffle_yeti_config.sh\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}".format(filepath))


def write_shuffling_script(experiment, data_file, save_dir, save_name):
    write_shuffled_data_generating_script(experiment, data_file, save_dir, save_name)
    write_shuffling_yeti_script(experiment)
    write_shuffling_submit_script(experiment)
    logger.info("done writing {} yeti scripts\n".format(experiment))


def setup_shuffle_model(conditions):
    for name, paths in conditions.items():
        logger.info("Creating working directory: {}".format(paths['shuffle_experiment']))
        os.makedirs(os.path.expanduser(paths['shuffle_experiment']))

        os.makedirs(os.path.expanduser(paths['shuffle_save_dir']), exist_ok=True)

        write_shuffling_script(paths['shuffle_experiment'],
                               paths['data_file'],
                               paths['shuffle_save_dir'],
                               paths['shuffle_save_name'])

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(paths['shuffle_experiment'])
        logger.debug("changed into dir: {}".format(os.getcwd()))

        shell_command = ".{}shuffle_start_job.sh".format(os.sep)
        logger.debug("About to run {}".format(shell_command))
        process_results = subprocess.run(shell_command, shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Shuffled dataset creation job submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def create_shuffle_configs(conditions, best_params):
    for name, paths in conditions.items():
        fname = os.path.join(paths['shuffle_experiment'], "write_shuffle_configs_for_loopy.m")
        with open(fname, 'w') as f:
            f.write("create_shuffle_config_files( ...\n")
            f.write("    'datapath', '{}.mat', ...\n".format(
                os.path.join(paths['shuffle_save_dir'], paths['shuffle_save_name'])))
            f.write("    'experiment_name', '{}', ...\n".format(paths['shuffle_experiment']))
            f.write("    'email_for_notifications', '{}', ...\n".format(EMAIL))
            f.write("    'yeti_user', '{}', ...\n".format(USER))
            f.write("    'compute_true_logZ', false, ...\n")
            f.write("    'reweight_denominator', 'mean_degree', ...\n")
            # f.write("    's_lambda_splits', 1, ...\n")
            # f.write("    's_lambdas_per_split', 1, ...\n")
            f.write("    's_lambda_min', {}, ...\n".format(best_params[name]['s_lambda']))
            f.write("    's_lambda_max', {}, ...\n".format(best_params[name]['s_lambda']))
            # f.write("    'density_splits', 1, ...\n")
            # f.write("    'densities_per_split', 1, ...\n")
            f.write("    'density_min', {}, ...\n".format(best_params[name]['density']))
            f.write("    'density_max', {}, ...\n".format(best_params[name]['density']))
            # f.write("    'p_lambda_splits', 1, ...\n")
            # f.write("    'p_lambdas_per_split', 1, ...\n")
            f.write("    'p_lambda_min', {}, ...\n".format(best_params[name]['p_lambda']))
            f.write("    'p_lambda_max', {}, ...\n".format(best_params[name]['p_lambda']))
            f.write("    'time_span', {}, ...\n".format(best_params[name]['time_span']))
            f.write("    'num_shuffle', {});\n".format(NSHUFFLE))
        f.closed
        logger.info("done writing {}".format(fname))

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(paths['shuffle_experiment'])
        logger.debug("changed into dir: {}".format(os.getcwd()))
        crf_util.run_matlab_command("write_shuffle_configs_for_{},".format(MODEL_TYPE),
                                    add_path=SOURCE_DIR)
        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def exec_shuffle_model(shuffle_experiment, **kwargs):
    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(shuffle_experiment)
    logger.debug("changed into dir: {}".format(os.getcwd()))

    process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
    if process_results.returncode:
        logger.critical("\nAre you on the yeti cluster? Job submission failed.")
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    logger.info("Shuffle jobs submitted.")

    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def simple_test_shuffle_datasets(shuffle_save_dir, shuffle_save_name, **kwargs):
    filebase = os.path.join(os.path.expanduser(shuffle_save_dir), shuffle_save_name + "_")
    return crf_util.get_max_job_done(filebase) >= NSHUFFLE


def test_shuffle_CRFs(shuffle_experiment, **kwargs):
    filebase = os.path.join(shuffle_experiment, "results", "result")
    return crf_util.get_max_job_done(filebase) >= NSHUFFLE


def exec_merge_shuffle_CRFs(shuffle_experiment, **kwargs):
    results_path = os.path.join(shuffle_experiment, "results")
    crf_util.run_matlab_command("save_and_merge_shuffled_models('{}'); ".format(results_path),
                                add_path=SOURCE_DIR)
    logger.info("Shuffle models merged and saved.\n")


def main(conditions):
    if conditions:
        if len(conditions) > 1:
            raise ValueError("Multiple conditions not currently supported.")
        conditions = get_conditions_metadata(conditions)
        # Create bare-bones shuffle folder and create shuffled datasets
        setup_shuffle_model(conditions)
        # Get best params
        best_params = {}
        for name, cond in conditions.items():
            best_params[name] = get_best_parameters(cond['experiment'])
        logger.info("Parameters for all conditions collected.\n")
        # create shuffle configs with best params (write and run write_configs_for_loopy.m)
        create_shuffle_configs(conditions, best_params)
        # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
        for name, meta in conditions.items():
            meta['to_test'] = simple_test_shuffle_datasets
            meta['to_run'] = exec_shuffle_model
        crf_util.wait_and_run(conditions)
        # Wait for shuffle CRFs to be done and run merge and save_shuffle
        for meta in conditions.values():
            meta['to_test'] = test_shuffle_CRFs
            meta['to_run'] = exec_merge_shuffle_CRFs
        crf_util.wait_and_run(conditions)
        # Extract ensemble neuron IDs. Write to disk?
    else:
        raise TypeError("At least one condition name must be passed on the command line.")


if __name__ == '__main__':
    start_time = time.time()

    # Each condition is a dict containing condition specific filepaths
    conditions = {name: {} for name in sys.argv[1:]}
    main(conditions)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
