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
logger.setLevel(logging.INFO)
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"
# *** end constants ***


def get_conditions_metadata(conditions):
    parameters_parser = crf_util.get_raw_configparser()
    parameters = crf_util.get_GridsearchOptions(parser=parameters_parser)
    parameters.update(crf_util.get_OtherOptions(parser=parameters_parser))
    parameters.update(crf_util.get_section_options('GeneralOptions', parser=parameters_parser))
    parameters.update(crf_util.get_section_options('YetiOptions', parser=parameters_parser))
    parameters.update(crf_util.get_section_options('YetiGenerateShuffledOptions',
                                                   parser=parameters_parser))
    for name, cond in conditions.items():
        cond.update(parameters)
        experiment = "{}_{}_{}".format(cond['experiment_name'], name, MODEL_TYPE)
        metadata = {'data_file': "{}_{}.mat".format(cond['experiment_name'], name),
                    'experiment': experiment,
                    'shuffle_save_name': "shuffled_{}_{}".format(cond['experiment_name'], name),
                    'shuffle_save_dir': os.path.join(cond['data_directory'],
                                                     "shuffled", experiment),
                    'shuffle_experiment': "shuffled_{}".format(experiment)
                    }
        metadata['shuffle_experiment_dir'] = os.path.join(
            cond['source_directory'], "fwMatch-darpa", "expt", metadata['shuffle_experiment'])
        conditions[name].update(metadata)
    return conditions


def write_shuffled_data_generating_script(params):
    # script for generate shuffled data
    filepath = os.path.join(params['shuffle_experiment'], "gn_shuff_data.m")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("addpath(genpath('{}'));\n".format(params['source_directory']))
        f.write("load(['{}']);\n".format(os.path.join(params['data_directory'],
                                                      params['data_file'])))
        f.write("fprintf('Loaded: %s\\n', ['{}']);\n".format(
            os.path.join(params['data_directory'], params['data_file'])))
        f.write("if exist('stimuli', 'var') ~= 1\n")
        f.write("    stimuli = [];\n")
        f.write("end\n")
        f.write("data_raw = data';\n")
        f.write("for i = 1:{}\n".format(params['num_shuffle']))
        f.write("\tdata = shuffle(data_raw,'exchange')';\n")
        f.write("\tsave(['{}_' num2str(i) '.mat'],'data','stimuli');\n".format(
            os.path.join(params['shuffle_save_dir'], params['shuffle_save_name'])))
        f.write("end\n")
        f.write("fprintf('done shuffling data\\n');\n")
    f.closed
    logger.info("done writing {}".format(filepath))


def write_shuffling_yeti_script(params):
    # write yeti script
    filepath = os.path.join(params['shuffle_experiment'], "shuffle_yeti_config.sh")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#shuffle_yeti_config.sh\n")
        f.write("#PBS -N {}\n".format(params['shuffle_experiment']))
        f.write("#PBS -W group_list={}\n".format(params['group_id']))
        f.write("#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
            params['yeti_gen_sh_nodes'], params['yeti_gen_sh_ppn'],
            params['yeti_gen_sh_walltime'], params['yeti_gen_sh_mem']))
        if params['email_notification'] == "num_jobs":
            # Always only 1 job, fully notify
            f.write("#PBS -m abe\n")
        else:
            # Use email_notification setting verbatim
            f.write("#PBS -m {}\n".format(params['email_notification']))
        f.write("#PBS -M {}\n".format(params['email']))
        f.write("#PBS -V\n")
        f.write("#set output and error directories (SSCC example here)\n")
        log_folder = "{}/".format(os.path.join(params['shuffle_experiment_dir'], "yeti_logs"))
        os.makedirs(os.path.expanduser(log_folder), exist_ok=True)
        f.write("#PBS -o localhost:{}\n".format(log_folder))
        f.write("#PBS -e localhost:{}\n".format(log_folder))
        f.write("#Command below is to execute Matlab code for Job Array (Example 4) " +
                "so that each part writes own output\n")
        f.write("matlab -nodesktop -nodisplay -r \"dbclear all;" +
                " addpath('{}');".format(os.path.join(params['shuffle_experiment_dir'])) +
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


def write_shuffling_script(params):
    write_shuffled_data_generating_script(params)
    write_shuffling_yeti_script(params)
    write_shuffling_submit_script(params['shuffle_experiment'])
    logger.info("done writing {} yeti scripts\n".format(params['shuffle_experiment']))


def setup_shuffle_model(conditions):
    for name, params in conditions.items():
        logger.info("Creating working directory: {}".format(params['shuffle_experiment']))
        os.makedirs(os.path.expanduser(params['shuffle_experiment']))

        os.makedirs(os.path.expanduser(params['shuffle_save_dir']), exist_ok=True)

        write_shuffling_script(params)

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(params['shuffle_experiment'])
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
    for name, params in conditions.items():
        fname = os.path.join(params['shuffle_experiment'], "write_shuffle_configs_for_loopy.m")
        with open(fname, 'w') as f:
            f.write("create_shuffle_config_files( ...\n")
            f.write("    'datapath', '{}.mat', ...\n".format(
                os.path.join(params['shuffle_save_dir'], params['shuffle_save_name'])))
            f.write("    'experiment_name', '{}', ...\n".format(params['shuffle_experiment']))
            f.write("    'email_for_notifications', '{}', ...\n".format(params['email']))
            f.write("    'yeti_user', '{}', ...\n".format(params['username']))
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
            f.write("    'num_shuffle', {});\n".format(params['num_shuffle']))
        f.closed
        logger.info("done writing {}".format(fname))

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(params['shuffle_experiment'])
        logger.debug("changed into dir: {}".format(os.getcwd()))
        crf_util.run_matlab_command("write_shuffle_configs_for_{},".format(MODEL_TYPE),
                                    add_path=params['source_directory'])
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


def simple_test_shuffle_datasets(shuffle_save_dir, shuffle_save_name, num_shuffle, **kwargs):
    filebase = os.path.join(os.path.expanduser(shuffle_save_dir), shuffle_save_name + "_")
    return crf_util.get_max_job_done(filebase) >= num_shuffle


def test_shuffle_CRFs(shuffle_experiment, num_shuffle, **kwargs):
    filebase = os.path.join(shuffle_experiment, "results", "result")
    return crf_util.get_max_job_done(filebase) >= num_shuffle


def exec_merge_shuffle_CRFs(shuffle_experiment, source_directory, **kwargs):
    results_path = os.path.join(shuffle_experiment, "results")
    crf_util.run_matlab_command("save_and_merge_shuffled_models('{}'); ".format(results_path),
                                add_path=source_directory)
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
