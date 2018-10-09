#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import crf_util
import gridsearch_train_crfs
import shuffled_control_crfs

import logging

logger = logging.getLogger("top")
logger.setLevel(logging.DEBUG)
init_log_handler = logging.StreamHandler(stream=sys.stdout)
init_log_handler.setLevel(logging.INFO)
logger.addHandler(init_log_handler)


def setup_logging(verbosity, debug_filelogging, expt_dir, **_):
    """Set up logging to monitor progress.

    All model module loggers only emit messages to their specific log files, and then pass the
    records up to the "top" logger. If you want messages to appear on screen, implementing a
    StreamHandler for the "top" logger here is necessary.

    Args:
        verbosity (int): How many messages to print to screen.
        debug_filelogging (bool): If True, print maximal messages to the log file.
        expt_dir (str (path)): Path at which to create the log file.
        **_: Catches and ignores any other keys in the passed dict, for convenience.
    """
    logger.addHandler(crf_util.get_StreamHandler(verbosity))
    logger.removeHandler(init_log_handler)
    logger.debug("Logging stream handler to sys.stdout added.")

    # Set log file to script's filename with a .log extension, in expt folder
    log_fname = os.path.join(os.path.expanduser(expt_dir), __file__) + ".log"
    logger.addHandler(
        crf_util.get_FileHandler(log_fname, debug_filelogging=debug_filelogging)
    )
    logger.debug("Logging file handler to {} added.".format(log_fname))


def merge_and_get_parameters(experiment, **kwargs):
    gridsearch_train_crfs.merge_save_train_models(experiment, kwargs["source_directory"])
    return gridsearch_train_crfs.get_best_parameters(experiment)


if __name__ == "__main__":
    start_time = time.time()

    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    params = {}
    if len(sys.argv) > 2:
        ini_fname = sys.argv[2]
        logger.debug("Received custom settings filename: {}".format(ini_fname))
        params.update(gridsearch_train_crfs.get_conditions_metadata(condition, ini_fname))
        params.update(shuffled_control_crfs.get_conditions_metadata(condition, ini_fname))
    else:
        params.update(gridsearch_train_crfs.get_conditions_metadata(condition))
        params.update(shuffled_control_crfs.get_conditions_metadata(condition))
    setup_logging(**params)

    gridsearch_train_crfs.setup_exec_train_model(params)
    # Create bare-bones shuffle folder and create shuffled datasets
    shuffled_control_crfs.setup_shuffle_model(params)
    # Wait for train CRF to be done
    # Run merge and save_best, grabbing best params
    params["to_test"] = params["test_gs_get_best_params"]
    params["to_run"] = merge_and_get_parameters
    best_params = crf_util.wait_and_run(params)
    logger.info("Parameters for all conditions collected.\n")
    # create shuffle configs with best params (write and run write_configs_for_loopy.m)
    shuffled_control_crfs.create_shuffle_configs(params, best_params)
    # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
    params["to_test"] = shuffled_control_crfs.simple_test_shuffle_datasets
    params["to_run"] = params["train_controls"]
    crf_util.wait_and_run(params)
    # Wait for shuffle CRFs to be done and run merge and save_shuffle
    params["to_test"] = shuffled_control_crfs.test_shuffle_CRFs
    params["to_run"] = shuffled_control_crfs.exec_merge_shuffle_CRFs
    crf_util.wait_and_run(params)
    # Extract ensemble neuron IDs. Write to disk?

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
