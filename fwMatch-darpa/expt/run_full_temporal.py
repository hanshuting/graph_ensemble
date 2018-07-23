#!/usr/bin/python
# -*- coding: utf-8 -*-
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
TRAIN_TEMPLATE_FOLDER_NAME = "{}_template".format(EXPT_NAME)
SHUFFLE_TEMPLATE_FOLDER_NAME = "shuffled_{}".format(TRAIN_TEMPLATE_FOLDER_NAME)
# *** end constants ***


def check_templates():
    """Make template folders from master templates if they do not exist already.
    """
    # Check train
    if not os.path.exists(TRAIN_TEMPLATE_FOLDER_NAME):
        shutil.copytree("temporal_template", TRAIN_TEMPLATE_FOLDER_NAME)

    if not os.path.exists(SHUFFLE_TEMPLATE_FOLDER_NAME):
        shutil.copytree("shuffled_template", SHUFFLE_TEMPLATE_FOLDER_NAME)


def setup_exec_train_model(condition_names):
    """Mostly follows old create_script.pl.

    Args:
        condition_names ([str]): List of condition names to setup.
    """
    for condition in condition_names:
        data_file = "{}_{}.mat".format(EXPT_NAME, condition)
        experiment = "{}_{}_{}".format(EXPT_NAME, condition, MODEL_TYPE)
        logger.info("Copying {} to {}".format(TRAIN_TEMPLATE_FOLDER_NAME, experiment))
        shutil.copytree(TRAIN_TEMPLATE_FOLDER_NAME, experiment)

        with open("{}{}write_configs_for_loopy.m".format(experiment, os.sep), 'w') as f:
            f.write("create_config_files( ...\n")
            f.write("    'datapath', '{}{}', ...\n".format(DATA_DIR, data_file))
            f.write("    'experiment_name', '{}', ...\n".format(experiment))
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
        logger.info("done writing write_configs_for_loopy.m")

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(experiment)
        logger.debug("changed into dir: {}".format(os.getcwd()))
        scommand = ("matlab -nodesktop -nodisplay -r \"try, write_configs_for_" +
                    "{}, catch, end, exit\"".format(MODEL_TYPE))
        logger.debug("About to run:\n{}".format(scommand))
        sargs = shlex.split(scommand)
        process_results = subprocess.run(sargs)
        if process_results.returncode:
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("\nTraining configs generated.")

        process_results = subprocess.run("./start_jobs.sh", shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Training job(s) submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def write_shuffled_data_generating_script(experiment, data_file, save_dir, save_name):
    # script for generate shuffled data
    filepath = "{}/gn_shuff_data.m".format(experiment)
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
        f.write("\tsave(['{}/{}_' num2str(i) '.mat'],'data','stimuli');\n".format(save_dir,
                                                                                  save_name))
        f.write("end\n")
        f.write("fprintf('done shuffling data\\n');\n")
    f.closed
    logger.info("done writing {}\n".format(filepath))


def write_shuffling_yeti_script(experiment):
    # write yeti script
    filepath = "{}/shuffle_yeti_config.sh".format(experiment)
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
                "{}fwMatch-darpa/expt/{}');gn_shuff_data; exit\"\n".format(SOURCE_DIR, experiment))
        f.write("#End of script\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}\n".format(filepath))


def write_shuffling_submit_script(experiment):
    # write submit script
    filepath = "{}/shuffle_start_job.sh".format(experiment)
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("qsub shuffle_yeti_config.sh\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}\n".format(filepath))


def write_shuffling_script(experiment, data_file, save_dir, save_name):
    write_shuffled_data_generating_script(experiment, data_file, save_dir, save_name)
    write_shuffling_yeti_script(experiment)
    write_shuffling_submit_script(experiment)
    logger.info("done writing {} yeti scripts\n".format(experiment))


def setup_shuffle_model(condition_names):
    for condition in condition_names:
        experiment = "shuffled_{}_{}_{}".format(EXPT_NAME, condition, MODEL_TYPE)
        logger.info("Copying {} to {}".format(SHUFFLE_TEMPLATE_FOLDER_NAME, experiment))
        shutil.copytree(SHUFFLE_TEMPLATE_FOLDER_NAME, experiment)
        data_file = "{}_{}.mat".format(EXPT_NAME, condition)

        save_dir = "{}shuffled/{}".format(DATA_DIR, experiment)
        # prev_umask = os.umask(mode=os.stat(DATA_DIR).st_mode)
        # os.makedirs(save_dir, mode=os.stat(DATA_DIR).st_mode, exist_ok=True)
        os.makedirs(os.path.expanduser(save_dir), exist_ok=True)
        save_name = "shuffled_{}_{}".format(EXPT_NAME, condition)

        write_shuffling_script(experiment, data_file, save_dir, save_name)

        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(experiment)
        logger.debug("changed into dir: {}".format(os.getcwd()))

        shell_command = "./shuffle_start_job.sh"
        logger.debug("About to run {}".format(shell_command))
        process_results = subprocess.run(shell_command, shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Shuffled dataset creation job submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def create_shuffle_configs(conditions):
    pass


if __name__ == '__main__':
    start_time = time.time()
    conditions = sys.argv[1:]
    if conditions:
        check_templates()
        setup_exec_train_model(conditions)
        # Create bare-bones shuffle folder and create shuffled datasets
        setup_shuffle_model(conditions)
        # Wait for train CRF to be done
        # Run merge and save_best, grabbing best params
        # create shuffle configs with best params (write and run write_configs_for_loopy.m)
        create_shuffle_configs(conditions)
        # Run shuffle/start_jobs.sh
        # Wait for shuffle CRFs to be done
        # Run merge and save_shuffle
        # Extract ensemble neuron IDs. Write to disk?
    else:
        raise TypeError("At least one condition name must be passed on the command line.")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
