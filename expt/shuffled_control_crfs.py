#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import subprocess

import crf_util
import gridsearch_train_crfs
import yeti_support

import logging
logger = logging.getLogger("top." + __name__)
logger.setLevel(logging.DEBUG)

# *** start constants ***
MODEL_TYPE = "loopy"
# *** end constants ***


def start_logfile(debug_filelogging, shuffle_experiment_dir, **_):
    log_fname = os.path.join(os.path.expanduser(shuffle_experiment_dir),
                             "shuffled_control_crfs.log")
    logfile_handler = crf_util.get_FileHandler(log_fname, debug_filelogging=debug_filelogging)
    logger.addHandler(logfile_handler)
    logger.debug("Logging file handler to {} added.".format(log_fname))


def get_conditions_metadata(condition):
    """Reads in settings.

    Args:
        condition (str): Condition name.

    Returns:
        dict: Metadata for condition.
    """
    parameters_parser = crf_util.get_raw_configparser()
    params = crf_util.get_GeneralOptions(parser=parameters_parser)
    experiment = "{}_{}_{}".format(params['experiment_name'], condition, MODEL_TYPE)
    metadata = {'data_file': "{}_{}.mat".format(params['experiment_name'], condition),
                'experiment': experiment,
                'shuffle_save_name': "shuffled_{}_{}".format(params['experiment_name'], condition),
                'shuffle_save_dir': os.path.join(params['data_directory'],
                                                 "shuffled", experiment),
                'shuffle_experiment': "shuffled_{}".format(experiment)
                }
    metadata['shuffle_experiment_dir'] = os.path.join(
        params['source_directory'], "expt", metadata['shuffle_experiment'])
    params.update(metadata)
    params["create_shuffles"] = create_shuffles
    # There is no training prep by default, so we set no-op function as placeholder
    params["shuffle_training_prep"] = lambda params: None
    params["train_controls"] = exec_shuffle_model

    # Update settings for cluster specified, if any
    if params["cluster_architecture"] == "yeti":
        logger.info("Yeti cluster architecture selected for shuffled controls.")
        params.update(yeti_support.get_yeti_shuff_metadata())

    return params


def create_shuffles(params):
    matlab_cmd = "gn_shuff_data('{}', '{}', {});".format(
        os.path.join(params['data_directory'], params['data_file']),
        params['shuffle_save_dir'],
        params['num_shuffle']
    )
    crf_util.run_matlab_command(matlab_cmd, add_path=params["source_directory"])


def setup_shuffle_model(params):
    logger.info("Creating working directory: {}".format(params['shuffle_experiment']))
    os.makedirs(os.path.expanduser(params['shuffle_experiment']))
    start_logfile(**params)

    os.makedirs(os.path.expanduser(params['shuffle_save_dir']), exist_ok=True)
    # TODO: Either clear pre-existing shuffled datasets, or skip regenerating any already there

    params["create_shuffles"](params)


def create_shuffle_configs(params, best_params):
    fname = os.path.join(params['shuffle_experiment'], "write_shuffle_configs_for_loopy.m")
    with open(fname, 'w') as f:
        f.write("create_shuffle_configs( ...\n")
        f.write("    'datapath', '{}.mat', ...\n".format(
            os.path.join(params['shuffle_save_dir'], params['shuffle_save_name'])))
        f.write("    'experiment_name', '{}', ...\n".format(params['shuffle_experiment']))
        try:
            f.write("    'email_for_notifications', '{}', ...\n".format(params['email']))
        except KeyError:
            logger.debug("No notifications email setting provided. Skipping.")
        try:
            f.write("    'yeti_user', '{}', ...\n".format(params['username']))
        except KeyError:
            logger.debug("No yeti username provided. Skipping.")
        f.write("    'compute_true_logZ', false, ...\n")
        f.write("    'reweight_denominator', 'mean_degree', ...\n")
        # f.write("    's_lambda_splits', 1, ...\n")
        # f.write("    's_lambdas_per_split', 1, ...\n")
        f.write("    's_lambda_min', {}, ...\n".format(best_params['s_lambda']))
        f.write("    's_lambda_max', {}, ...\n".format(best_params['s_lambda']))
        # f.write("    'density_splits', 1, ...\n")
        # f.write("    'densities_per_split', 1, ...\n")
        f.write("    'density_min', {}, ...\n".format(best_params['density']))
        f.write("    'density_max', {}, ...\n".format(best_params['density']))
        # f.write("    'p_lambda_splits', 1, ...\n")
        # f.write("    'p_lambdas_per_split', 1, ...\n")
        f.write("    'p_lambda_min', {}, ...\n".format(best_params['p_lambda']))
        f.write("    'p_lambda_max', {}, ...\n".format(best_params['p_lambda']))
        f.write("    'edges', '{}', ...\n".format(params['edges'].lower()))
        f.write("    'no_same_neuron_edges', {}, ...\n".format(
            str(params['no_same_neuron_edges']).lower()))
        f.write("    'time_span', {}, ...\n".format(best_params['time_span']))
        f.write("    'num_shuffle', {});\n".format(params['num_shuffle']))
    f.closed
    logger.info("done writing {}".format(fname))

    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(params['shuffle_experiment'])
    logger.debug("changed into dir: {}".format(os.getcwd()))
    crf_util.run_matlab_command("write_shuffle_configs_for_{},".format(MODEL_TYPE),
                                add_path=params['source_directory'])

    params["shuffle_training_prep"](params)

    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def exec_shuffle_model(source_directory, num_shuffle, **kwargs):
    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(source_directory)
    logger.debug("changed into dir: {}".format(os.getcwd()))
    for cur_shuffle in range(1, num_shuffle + 1):
        process_results = subprocess.run(".{}run.sh {}".format(os.sep, cur_shuffle), shell=True)
        if process_results.returncode:
            logger.critical("Training failed on shuffled control {}".format(cur_shuffle))
            raise RuntimeError("Received non-zero return code: {}".format(process_results))
        logger.info("Trained shuffled control model {} out of {}.".format(cur_shuffle,
                                                                          num_shuffle))
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


def main(condition):
    """Summary

    Args:
        condition (str): Condition name.
    """
    params = get_conditions_metadata(condition)

    # Update logging if module is invoked from the command line
    if __name__ == '__main__':
        # Assume top log position
        logger = logging.getLogger("top")
        logger.setLevel(logging.DEBUG)
        # Create stdout log handler
        verbosity = params['verbosity']
        logger.addHandler(crf_util.get_StreamHandler(verbosity))
        logger.debug("Logging stream handler to sys.stdout added.")

    # Create bare-bones shuffle folder and create shuffled datasets
    setup_shuffle_model(params)
    # Get best params
    best_params = gridsearch_train_crfs.get_best_parameters(params['experiment'])
    logger.info("Parameters for {} collected.\n".format(condition))
    # create shuffle configs with best params (write and run write_configs_for_loopy.m)
    create_shuffle_configs(params, best_params)
    # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
    params['to_test'] = simple_test_shuffle_datasets
    params['to_run'] = params["train_controls"]
    crf_util.wait_and_run(params)
    # Wait for shuffle CRFs to be done and run merge and save_shuffle
    params['to_test'] = test_shuffle_CRFs
    params['to_run'] = exec_merge_shuffle_CRFs
    crf_util.wait_and_run(params)
    # Extract ensemble neuron IDs. Write to disk?


if __name__ == '__main__':
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    main(condition)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
