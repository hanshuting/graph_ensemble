#!/usr/bin/python
# -*- coding: utf-8 -*-
import time
import sys
import os
import shutil
import shlex
import subprocess

# *** START USER EDITABLE VARIABLES ***
EXPT_NAME = "temporal"
DATA_DIR = "~/data/"
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
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"
TRAIN_TEMPLATE_FOLDER_NAME = "{}_template".format(EXPT_NAME)
#SHUFFLE_TEMPLATE_FOLDER_NAME = # TODO
# *** end constants ***


def check_templates():
    """Make template folders from master templates if they do not exist already.
    """
    # Check train
    if not os.path.exists(TRAIN_TEMPLATE_FOLDER_NAME):
        shutil.copytree("temporal_template", TRAIN_TEMPLATE_FOLDER_NAME)

    # TODO: Check stuffle


def setup_exec_train_model(condition_names):
    """Mostly follows old create_script.pl.

    Args:
        condition_names ([str]): List of condition names to setup.
    """
    for condition in condition_names:
        data_file = "{}_{}.mat".format(EXPT_NAME, condition)
        experiment = "{}_{}_{}".format(EXPT_NAME, condition, MODEL_TYPE)
        print("Copying {} to {}".format(TRAIN_TEMPLATE_FOLDER_NAME, experiment))
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
        print("done writing write_configs_for_loopy.m")

        curr_dir = os.getcwd()
        print("curr_dir = {}.".format(curr_dir))
        os.chdir(experiment)
        print("changed into dir: {}".format(os.getcwd()))
        scommand = ("matlab -nodesktop -nodisplay -r \"try, write_configs_for_" +
                    "{}, catch, end, exit\"".format(MODEL_TYPE))
        print("About to run:\n{}".format(scommand))
        sargs = shlex.split(scommand)
        process_results = subprocess.run(sargs)
        if process_results.returncode:
            raise RuntimeError("Received non-zero return code: {}".format(process_results))

        process_results = subprocess.run("./start_jobs.sh", shell=True)
        if process_results.returncode:
            print("Are you on the yeti cluster?")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))

        os.chdir(curr_dir)
        print("changed back to dir: {}".format(os.getcwd()))


if __name__ == '__main__':
    start_time = time.time()
    conditions = sys.argv[1:]
    if conditions:
        check_templates()
        setup_exec_train_model(conditions)
        # Create bare-bones shuffle folder
        # run shuffle dataset creation, if needed
        # Wait for train CRF to be done
        # Run merge and save_best, grabbing best params
        # create shuffle configs with best params
        # Run shuffle/start_jobs.sh
        # Wait for shuffle CRFs to be done
        # Run merge and save_shuffle
        # Extract ensemble neuron IDs. Write to disk?
    else:
        raise TypeError("At least one condition name must be passed on the command line.")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
