#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import subprocess
import crf_util

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))

# *** START USER EDITABLE VARIABLES ***
logger.setLevel(logging.INFO)
EXPT_NAME = "experiment"
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
EXPT_DIR = os.path.join(SOURCE_DIR, "fwMatch-darpa", "expt")

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


def setup_exec_train_model(conditions):
    """Mostly follows old create_script.pl.

    Args:
        conditions (dict): Dict of condition names: paths.
    """
    for name, paths in conditions.items():
        logger.info("Creating working directory: {}".format(paths['experiment']))
        os.makedirs(os.path.expanduser(paths['experiment']))

        fname = os.path.join(paths['experiment'], "write_configs_for_loopy.m")
        # TODO: Just call create_config_files directly
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
        crf_util.run_matlab_command("try, write_configs_for_{}, catch, end,".format(MODEL_TYPE),
                                    add_path=SOURCE_DIR)
        logger.info("\nTraining configs generated.")

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Training job(s) submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def merge_save_train_models(experiment, **kwargs):
    results_path = os.path.join(experiment, "results")
    return crf_util.run_matlab_command("save_best_params('{}'); ".format(results_path),
                                       add_path=SOURCE_DIR)


def get_best_parameters(experiment, **kwargs):
    """Get the best parameters from gridsearch results.

    Args:
        experiment (str): Path to working directory containing results subdirectory.

    Returns:
        dict: Best parameters from the gridsearch. Parameters stored as a dict with
            PARAMS_TO_EXTRACT as the keys.
    """
    best_params = {}

    # grab and return best params
    results_path = os.path.join(experiment, "results")
    with open(os.path.join(results_path, "best_parameters.txt"), 'r') as f:
            for param in PARAMS_TO_EXTRACT:
                best_params[param] = float(f.readline())
    f.closed

    return best_params


def test_train_CRFs(experiment, **kwargs):
    filebase = os.path.join(experiment, "results", "result")
    num_jobs = 1
    for param in [S_LAMBDAS, DENSITIES, P_LAMBDAS]:
        num_jobs *= param['num_points'] if param['parallize'] else 1
    return crf_util.get_max_job_done(filebase) >= num_jobs


def main(conditions):
    """Summary

    Args:
        conditions (dict): Each key refers to a condition, and its value is a dict containing
            condition specific filepaths.
    """
    if conditions:
        if len(conditions) > 1:
            raise ValueError("Multiple conditions not currently supported.")
        conditions = get_conditions_metadata(conditions)
        setup_exec_train_model(conditions)
        # Wait for train CRF to be done
        # Run merge and save_best
        for cond in conditions.values():
            cond['to_test'] = test_train_CRFs
            cond['to_run'] = merge_save_train_models
        crf_util.wait_and_run(conditions)
        for cond in conditions.values():
            best_params_path = os.path.join(EXPT_DIR, cond['experiment'],
                                            "results", "best_parameters.txt")
            logger.info("Grid search complete. Best parameters in {}".format(best_params_path) +
                        " in the following order:\n{}\n".format(PARAMS_TO_EXTRACT))
    else:
        raise TypeError("At least one condition name must be passed on the command line.")


if __name__ == '__main__':
    start_time = time.time()

    conditions = {name: {} for name in sys.argv[1:]}
    main(conditions)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
