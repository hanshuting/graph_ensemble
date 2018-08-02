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

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))

# *** START USER EDITABLE VARIABLES ***
logger.setLevel(logging.INFO)
# *** END USER EDITABLE VARIABLES ***

# *** start constants ***
MODEL_TYPE = "loopy"

# These parameters and their order must match best_parameters.txt.
# See save_best_parameters.m for best_parameters.txt creation.
PARAMS_TO_EXTRACT = ['s_lambda', 'density', 'p_lambda', 'time_span']
# *** end constants ***


def get_conditions_metadata(conditions):
    parameters_parser = crf_util.get_raw_configparser()
    parameters = crf_util.get_GridsearchOptions(parser=parameters_parser)
    parameters.update(crf_util.get_OtherOptions(parser=parameters_parser))
    parameters.update(crf_util.get_section_options('GeneralOptions', parser=parameters_parser))
    parameters.update(crf_util.get_section_options('YetiOptions', parser=parameters_parser))
    for name, cond in conditions.items():
        cond.update(parameters)
        metadata = {'data_file': "{}_{}.mat".format(cond['experiment_name'], name),
                    'experiment': "{}_{}_{}".format(cond['experiment_name'], name, MODEL_TYPE),
                    'expt_dir': os.path.join(cond['source_directory'], "fwMatch-darpa", "expt")
                    }
        cond.update(metadata)
    return conditions


def create_write_configs_for_loopy_m(name, params):
        logger.info("Creating working directory: {}".format(params['experiment']))
        os.makedirs(os.path.expanduser(params['experiment']))

        fname = os.path.join(params['experiment'], "write_configs_for_loopy.m")
        # TODO: Just call create_config_files directly
        with open(fname, 'w') as f:
            f.write("create_config_files( ...\n")
            f.write("    'datapath', '{}{}', ...\n".format(params['data_directory'],
                                                           params['data_file']))
            f.write("    'experiment_name', '{}', ...\n".format(params['experiment']))
            f.write("    'email_for_notifications', '{}', ...\n".format(params['email']))
            f.write("    'yeti_user', '{}', ...\n".format(params['username']))
            f.write("    'compute_true_logZ', false, ...\n")
            f.write("    'reweight_denominator', 'mean_degree', ...\n")
            s_lambdas = params['S_LAMBDAS']
            f.write("    's_lambda_splits', {}, ...\n".format(
                s_lambdas['num_points'] if s_lambdas['parallize'] else 1))
            f.write("    's_lambdas_per_split', {}, ...\n".format(
                1 if s_lambdas['parallize'] else s_lambdas['num_points']))
            f.write("    's_lambda_min', {}, ...\n".format(s_lambdas['min']))
            f.write("    's_lambda_max', {}, ...\n".format(s_lambdas['max']))
            densities = params['DENSITIES']
            f.write("    'density_splits', {}, ...\n".format(
                densities['num_points'] if densities['parallize'] else 1))
            f.write("    'densities_per_split', {}, ...\n".format(
                1 if densities['parallize'] else densities['num_points']))
            f.write("    'density_min', {}, ...\n".format(densities['min']))
            f.write("    'density_max', {}, ...\n".format(densities['max']))
            p_lambdas = params['P_LAMBDAS']
            f.write("    'p_lambda_splits', {}, ...\n".format(
                p_lambdas['num_points'] if p_lambdas['parallize'] else 1))
            f.write("    'p_lambdas_per_split', {}, ...\n".format(
                1 if p_lambdas['parallize'] else p_lambdas['num_points']))
            f.write("    'p_lambda_min', {}, ...\n".format(p_lambdas['min']))
            f.write("    'p_lambda_max', {}, ...\n".format(p_lambdas['max']))
            f.write("    'time_span', {});\n".format(params['time_span']))
    f.closed
    logger.info("done writing {}".format(fname))


def create_start_jobs_sh(experiment):
    # Expect to be in the experiment folder already when writing this
    # TODO: yeti specific
    fname = "start_jobs.sh"
    with open(fname, 'w') as f:
        # Clear out any basic remainders from previous runs
        f.write("rm -f ./results/result*.mat\n")
        f.write("rm -f ./yeti_logs/*\n")
        f.write("rm -f ./job_logs/*\n")
        f.write("cd ../.. && qsub {}\n".format(os.path.join("expt", experiment, "yeti_config.sh")))
    f.closed

    # make sure file is executable
    os.chmod(fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode)
    logger.info("Created " + os.path.join(experiment, fname) + ".")


def setup_exec_train_model(conditions):
    """Mostly follows old create_script.pl.

        Expect to be in expt folder at start.

    Args:
        conditions (dict): Dict of condition names: params.
    """
    for name, params in conditions.items():
        create_write_configs_for_loopy_m(name, params)

        # Move into experiment folder
        curr_dir = os.getcwd()
        logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(params['experiment'])
        logger.debug("changed into dir: {}".format(os.getcwd()))
        crf_util.run_matlab_command("try, write_configs_for_{}, catch, end,".format(MODEL_TYPE),
                                    add_path=params['source_directory'])
        logger.info("\nTraining configs generated.")

        create_start_jobs_sh(params['experiment_name'])

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Training job(s) submitted.")

        os.chdir(curr_dir)
        logger.debug("changed back to dir: {}".format(os.getcwd()))


def merge_save_train_models(experiment, source_directory, **kwargs):
    results_path = os.path.join(experiment, "results")
    return crf_util.run_matlab_command("save_best_params('{}'); ".format(results_path),
                                       add_path=source_directory)


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


def test_train_CRFs(experiment, S_LAMBDAS, DENSITIES, P_LAMBDAS, **kwargs):
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
            best_params_path = os.path.join(cond['expt_dir'], cond['experiment'],
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
