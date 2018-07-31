#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
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
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"

# These parameters and their order must match best_parameters.txt.
# See save_best_parameters.m for best_parameters.txt creation.
PARAMS_TO_EXTRACT = ['s_lambda', 'density', 'p_lambda', 'time_span']
# *** end constants ***


def get_conditions_metadata(conditions):
    for name in conditions:
        metadata = {'data_file': "{}_{}.mat".format(EXPT_NAME, name),
                    'experiment': "{}_{}_{}".format(EXPT_NAME, name, MODEL_TYPE),
                    }
        conditions[name].update(metadata)
    return conditions


def run_command(scommand):
    logger.debug("About to run:\n{}".format(scommand))
    sargs = shlex.split(scommand)
    process_results = subprocess.run(sargs)
    if process_results.returncode:
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    return process_results


def run_matlab_command(scommand):
    return run_command("matlab -nodesktop -nodisplay -nosplash -r \"" +
                       "addpath(genpath('{}')); ".format(SOURCE_DIR) +
                       scommand +
                       "exit\"")


def setup_exec_train_model(conditions):
    """Mostly follows old create_script.pl.

    Args:
        conditions (dict): Dict of condition names: paths.
    """
    for name, paths in conditions.items():
        logger.info("Creating working directory: {}".format(paths['experiment']))
        os.makedirs(os.path.expanduser(paths['experiment']))

        fname = "{}{}write_configs_for_loopy.m".format(paths['experiment'], os.sep)
        # TODO: Just call write_configs_for_loopy directly
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
        run_command("matlab -nodesktop -nodisplay -r \"" +
                    "addpath(genpath('{}')); ".format(SOURCE_DIR) +
                    "try, write_configs_for_{}, catch, end, exit\"".format(MODEL_TYPE))
        logger.info("\nTraining configs generated.")

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Training job(s) submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def merge_save_train_models(experiment, **kwargs):
    results_path = "{0}{1}results{1}".format(experiment, os.sep)
    return run_matlab_command("save_best_params('{}'); ".format(results_path))


def get_best_parameters(experiment, **kwargs):
    """Extract the best parameters from gridsearch results.

    Args:
        experiment (str): Path to working directory containing results subdirectory.

    Returns:
        dict: Best parameters from the gridsearch. Parameters stored as a dict with
            PARAMS_TO_EXTRACT as the keys.
    """
    # merge & save models
    merge_save_train_models(experiment)

    best_params = {}
    results_path = "{0}{1}results{1}".format(experiment, os.sep)

    # grab and return best params
    with open(results_path + "best_parameters.txt", 'r') as f:
            for param in PARAMS_TO_EXTRACT:
                best_params[param] = float(f.readline())
    f.closed

    return best_params


def get_max_job_done(filebase, filesuffix=".mat"):
    filebase = os.path.expanduser(filebase)
    job = 1
    while os.path.exists("{}{}{}".format(filebase, job, filesuffix)):
        job += 1
    return job - 1


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


if __name__ == '__main__':
    start_time = time.time()

    # Each condition is a dict containing condition specific filepaths
    conditions = {name: {} for name in sys.argv[1:]}
    if conditions:
        if len(conditions) > 1:
            raise ValueError("Multiple conditions not currently supported.")
        conditions = get_conditions_metadata(conditions)
        setup_exec_train_model(conditions)
        # Wait for train CRF to be done
        # Run merge and save_best, grabbing best params
        for cond in conditions.values():
            cond['to_test'] = test_train_CRFs
            cond['to_run'] = merge_save_train_models
        wait_and_run(conditions)
        full_results_path = "{0}{1}{2}{1}results{1}".format(SOURCE_DIR, os.sep,
                                                            conditions['experiment'])
        logger.info("Grid search complete. Best parameters in {}".format(full_results_path) +
                    "best_parameters.txt in the following order:\n{}\n".format())
    else:
        raise TypeError("At least one condition name must be passed on the command line.")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
