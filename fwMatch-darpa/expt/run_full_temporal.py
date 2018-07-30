#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import stat
import shutil
import shlex
import subprocess

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
S_LAMBDAS = {'parallize': True,      # If True, each point generates a distinct config and job
             'num_points': 3,        # Number of parameter values to test between min and max
             'min': 2e-03,
             'max': 5e-01}
DENSITIES = {'parallize': False,      # If True, each point generates a distinct config and job
             'num_points': 3,        # Number of parameter values to test between min and max
             'min': 0.1,
             'max': 0.3}
P_LAMBDAS = {'parallize': True,      # If True, each point generates a distinct config and job
             'num_points': 2,        # Number of parameter values to test between min and max
             'min': 1e+01,
             'max': 1e+04}
TIME_SPAN = 2
NSHUFFLE = 100
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"

# These parameters and their order must match best_parameters.txt.
# See save_best_parameters.m for best_parameters.txt creation.
PARAMS_TO_EXTRACT = ['s_lambda', 'density', 'p_lambda', 'time_span']

TRAIN_TEMPLATE_FOLDER_NAME = "{}_template".format(EXPT_NAME)
SHUFFLE_TEMPLATE_FOLDER_NAME = "shuffled_{}".format(TRAIN_TEMPLATE_FOLDER_NAME)
# *** end constants ***


def get_conditions_metadata(conditions):
    for name in conditions:
        metadata = {'data_file': "{}_{}.mat".format(EXPT_NAME, name),
                    'experiment': "{}_{}_{}".format(EXPT_NAME, name, MODEL_TYPE),
                    'shuffle_save_name': "shuffled_{}_{}".format(EXPT_NAME, name)
                    }
        metadata['shuffle_save_dir'] = "{}shuffled{}{}".format(DATA_DIR, os.sep,
                                                               metadata['experiment'])
        metadata['shuffle_experiment'] = "shuffled_{}".format(metadata['experiment'])
        conditions[name].update(metadata)
    return conditions


def check_templates():
    """Make template folders from master templates if they do not exist already.
    """
    # Check train
    if not os.path.exists(TRAIN_TEMPLATE_FOLDER_NAME):
        shutil.copytree("temporal_template", TRAIN_TEMPLATE_FOLDER_NAME)

    if not os.path.exists(SHUFFLE_TEMPLATE_FOLDER_NAME):
        shutil.copytree("shuffled_template", SHUFFLE_TEMPLATE_FOLDER_NAME)


def run_command(scommand):
    logger.debug("About to run:\n{}".format(scommand))
    sargs = shlex.split(scommand)
    process_results = subprocess.run(sargs)
    if process_results.returncode:
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    return process_results


def setup_exec_train_model(conditions):
    """Mostly follows old create_script.pl.

    Args:
        conditions (dict): Dict of condition names: paths.
    """
    for name, paths in conditions.items():
        logger.info("Copying {} to {}".format(TRAIN_TEMPLATE_FOLDER_NAME, paths['experiment']))
        shutil.copytree(TRAIN_TEMPLATE_FOLDER_NAME, paths['experiment'])

        fname = "{}{}write_configs_for_loopy.m".format(paths['experiment'], os.sep)
        with open(fname, 'w') as f:
            f.write("create_config_files( ...\n")
            f.write("    'datapath', '{}{}', ...\n".format(DATA_DIR, paths['data_file']))
            f.write("    'experiment_name', '{}', ...\n".format(paths['experiment']))
            f.write("    'email_for_notifications', '{}', ...\n".format(EMAIL))
            f.write("    'yeti_user', '{}', ...\n".format(USER))
            f.write("    'compute_true_logZ', false, ...\n")
            f.write("    'reweight_denominator', 'mean_degree', ...\n")
            f.write("    's_lambda_splits', {}, ...\n".format(
                S_LAMBDAS['num_points'] if S_LAMBDAS['parallize'] else 1))
            f.write("    's_lambdas_per_split', {}, ...\n".format(
                1 if S_LAMBDAS['parallize'] else S_LAMBDAS['num_points']))
            f.write("    's_lambda_min', {}, ...\n".format(S_LAMBDAS['min']))
            f.write("    's_lambda_max', {}, ...\n".format(S_LAMBDAS['max']))
            f.write("    'density_splits', {}, ...\n".format(
                DENSITIES['num_points'] if DENSITIES['parallize'] else 1))
            f.write("    'densities_per_split', {}, ...\n".format(
                1 if DENSITIES['parallize'] else DENSITIES['num_points']))
            f.write("    'density_min', {}, ...\n".format(DENSITIES['min']))
            f.write("    'density_max', {}, ...\n".format(DENSITIES['max']))
            f.write("    'p_lambda_splits', {}, ...\n".format(
                P_LAMBDAS['num_points'] if P_LAMBDAS['parallize'] else 1))
            f.write("    'p_lambdas_per_split', {}, ...\n".format(
                1 if P_LAMBDAS['parallize'] else P_LAMBDAS['num_points']))
            f.write("    'p_lambda_min', {}, ...\n".format(P_LAMBDAS['min']))
            f.write("    'p_lambda_max', {}, ...\n".format(P_LAMBDAS['max']))
            f.write("    'time_span', {});\n".format(TIME_SPAN))
        f.closed
        logger.info("done writing {}".format(fname))

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(paths['experiment'])
        logger.debug("changed into dir: {}".format(os.getcwd()))
        run_command("matlab -nodesktop -nodisplay -r \"try, write_configs_for_" +
                    "{}, catch, end, exit\"".format(MODEL_TYPE))
        logger.info("\nTraining configs generated.")

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Training job(s) submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def write_shuffled_data_generating_script(experiment, data_file, save_dir, save_name):
    # script for generate shuffled data
    filepath = "{}{}gn_shuff_data.m".format(experiment, os.sep)
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
        f.write("\tsave(['{}{}{}_' num2str(i) '.mat'],'data','stimuli');\n".format(
            save_dir, os.sep, save_name))
        f.write("end\n")
        f.write("fprintf('done shuffling data\\n');\n")
    f.closed
    logger.info("done writing {}".format(filepath))


def write_shuffling_yeti_script(experiment):
    # write yeti script
    filepath = "{}{}shuffle_yeti_config.sh".format(experiment, os.sep)
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
        f.write("matlab -nodesktop -nodisplay -r \"dbclear all; addpath('" +
                "{0}fwMatch-darpa{2}expt{2}{1}');".format(SOURCE_DIR, experiment, os.sep) +
                "gn_shuff_data; exit\"\n")
        f.write("#End of script\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}".format(filepath))


def write_shuffling_submit_script(experiment):
    # write submit script
    filepath = "{}{}shuffle_start_job.sh".format(experiment, os.sep)
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
        logger.info("Copying {} to {}".format(SHUFFLE_TEMPLATE_FOLDER_NAME,
                                              paths['shuffle_experiment']))
        shutil.copytree(SHUFFLE_TEMPLATE_FOLDER_NAME, paths['shuffle_experiment'])

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
        fname = "{}{}write_configs_for_loopy.m".format(paths['shuffle_experiment'], os.sep)
        with open(fname, 'w') as f:
            f.write("create_shuffle_config_files( ...\n")
            f.write("    'datapath', '{}{}{}.mat', ...\n".format(paths['shuffle_save_dir'], os.sep,
                                                                 paths['shuffle_save_name']))
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
        run_command("matlab -nodesktop -nodisplay -r \"" +
                    "write_configs_for_{}, exit\"".format(MODEL_TYPE))

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def get_best_parameters(experiment, **kwargs):
    """Extract the best parameters from gridsearch results.

    Args:
        experiment (str): Path to working directory containing results subdirectory.

    Returns:
        dict: Best parameters from the gridsearch. Parameters stored as a dict with
            PARAMS_TO_EXTRACT as the keys.
    """
    best_params = {}
    results_path = "{0}{1}results{1}".format(experiment, os.sep)

    # merge & save models
    run_command("matlab -nodesktop -nodisplay -nosplash -r \"" +
                "addpath(genpath('{}')); ".format(SOURCE_DIR) +
                "save_best_params('{}'); ".format(results_path) +
                "exit\"")

    # grab and return best params
    with open(results_path + "best_parameters.txt", 'r') as f:
            for param in PARAMS_TO_EXTRACT:
                best_params[param] = float(f.readline())
    f.closed

    return best_params


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


def test_shuffle_datasets_exist(paths, **kwargs):
    job = 1
    while job <= NSHUFFLE:
        while os.path.exists("{}{}{}_{}.mat".format(paths['shuffle_save_dir'],
                                                    os.sep,
                                                    paths['shuffle_save_name'],
                                                    job)):
            job += 1
        yield job > NSHUFFLE
    return True


def simple_test_shuffle_datasets(shuffle_save_dir, shuffle_save_name, **kwargs):
    filebase = os.path.expanduser("{}{}{}".format(shuffle_save_dir, os.sep, shuffle_save_name))
    job = 1
    while os.path.exists("{}_{}.mat".format(filebase, job)):
        job += 1
    return job > NSHUFFLE


def get_max_job_done(filebase, filesuffix=".mat"):
    filebase = os.path.expanduser(filebase)
    job = 1
    while os.path.exists("{}{}{}".format(filebase, job, filesuffix)):
        job += 1
    return job - 1


def test_shuffle_CRFs(shuffle_experiment, **kwargs):
    filebase = "{0}{1}results{1}result".format(shuffle_experiment, os.sep)
    return get_max_job_done(filebase) >= NSHUFFLE


def test_train_CRFs(experiment, **kwargs):
    filebase = "{0}{1}results{1}result".format(experiment, os.sep)
    num_jobs = 1
    for param in [S_LAMBDAS, DENSITIES, P_LAMBDAS]:
        num_jobs *= param['num_points'] if param['parallize'] else 1
    return get_max_job_done(filebase) >= num_jobs


def wait_and_run(conditions_to_check, wait_seconds=5):
    """Execute specified functions after their corresponding tests pass, pausing between tests.

    Args:
        conditions_to_check (dict of dicts): An item per waiting task and subsequent execution.
            Expects each top-level key to have an associated dict with at least:
                'to_test': a function that returns true when testing should conclude and execution
                    should begin.
                'to_run': the function to run once 'to_test' returns true. Ret
            The full dict of each top level key is passed as kwargs to its 'to_test' and 'to_run'.
        wait_seconds (float, optional): Number of seconds to wait per 'to_test' iterations.
    """
    return_vals = {}
    conditions_remaining = {name: None for name in conditions_to_check}
    logger.debug("Start waiting for\n{}".format(conditions_to_check))
    num_waits = 0
    while conditions_remaining:
        stop_checking = []
        for name in conditions_remaining:
            to_check = conditions_to_check[name]
            if to_check['to_test'](**to_check):
                logger.debug("{}['to_test'] passed.".format(name))
                # TODO: Parallize here so we can run but still continue to test others?
                return_vals[name] = to_check['to_run'](**to_check)
                logger.info("{} for {} completed.".format(to_check['to_run'].__name__, name))
                logger.debug("return_vals['{}'] = {}".format(name, return_vals[name]))
                stop_checking.append(name)
        for finished in stop_checking:
            del conditions_remaining[finished]

        time.sleep(wait_seconds)
        num_waits += 1
        if (num_waits % 100) == 0:
            logger.info("Waited for {} sleep cycles so far. Currently waiting for:\n{}".format(
                num_waits, conditions_to_check))
        elif (num_waits % 20) == 0:
            logger.debug("Waited for {} sleep cycles so far. Currently waiting for:\n{}".format(
                num_waits, conditions_to_check))

    logger.debug("Done waiting for {}.\n".format(conditions_to_check.keys()))
    # TODO: Returning all values together means the last test to pass blocks returing others.
    return return_vals


def exec_merge_shuffle_CRFs(shuffle_experiment, **kwargs):
    results_path = "{0}{1}results{1}".format(shuffle_experiment, os.sep)
    run_command("matlab -nodesktop -nodisplay -nosplash -r \"" +
                "addpath(genpath('{}')); ".format(SOURCE_DIR) +
                "save_and_merge_shuffled_models('{}'); ".format(results_path) +
                "exit\"")
    logger.info("Shuffle models merged and saved.")


if __name__ == '__main__':
    start_time = time.time()

    # Each condition is a dict containing condition specific filepaths
    conditions = {name: {} for name in sys.argv[1:]}
    if conditions:
        if len(conditions) > 1:
            raise ValueError("Multiple conditions not currently supported.")
        conditions = get_conditions_metadata(conditions)
        check_templates()
        setup_exec_train_model(conditions)
        # Create bare-bones shuffle folder and create shuffled datasets
        setup_shuffle_model(conditions)
        # Wait for train CRF to be done
        # Run merge and save_best, grabbing best params
        for cond in conditions.values():
            cond['to_test'] = test_train_CRFs
            cond['to_run'] = get_best_parameters
        best_params = wait_and_run(conditions)
        logger.info("Parameters for all conditions collected.\n")
        # create shuffle configs with best params (write and run write_configs_for_loopy.m)
        create_shuffle_configs(conditions, best_params)
        # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
        for name, meta in conditions.items():
            meta['to_test'] = simple_test_shuffle_datasets
            meta['to_run'] = exec_shuffle_model
        wait_and_run(conditions)
        # Wait for shuffle CRFs to be done and run merge and save_shuffle
        for meta in conditions.values():
            meta['to_test'] = test_shuffle_CRFs
            meta['to_run'] = exec_merge_shuffle_CRFs
        wait_and_run(conditions)
        # Extract ensemble neuron IDs. Write to disk?
    else:
        raise TypeError("At least one condition name must be passed on the command line.")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
