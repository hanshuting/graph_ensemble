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
    stream_handler = logging.StreamHandler(stream=sys.stdout)
    stream_handler.setLevel(crf_util.loglevel_from_verbosity(verbosity))
    logger.addHandler(stream_handler)
    logger.debug("Logging stream handler to sys.stdout added.")

    # Set log file to script's filename with a .log extension, in expt folder
    log_fname = os.path.join(os.path.expanduser(expt_dir), __file__) + ".log"
    # Overwrite previous log file every time via mode='w'. Change to mode='a' to append.
    logfile_handler = logging.FileHandler(log_fname, mode='w')
    logfile_handler.setLevel(logging.DEBUG if debug_filelogging else logging.INFO)
    logfile_format = logging.Formatter('%(asctime)s - %(levelname)s@%(name)s - %(message)s')
    logfile_handler.setFormatter(logfile_format)
    logger.addHandler(logfile_handler)
    logger.debug("Logging file handler to {} added.".format(log_fname))


def merge_and_get_parameters(experiment, **kwargs):
    gridsearch_train_crfs.merge_save_train_models(experiment, kwargs['source_directory'])
    return gridsearch_train_crfs.get_best_parameters(experiment)


if __name__ == '__main__':
    start_time = time.time()

    # Each condition is a dict containing condition specific filepaths
    conditions = {name: {} for name in sys.argv[1:]}
    if conditions:
        if len(conditions) > 1:
            raise ValueError("Multiple conditions not currently supported.")
        conditions = gridsearch_train_crfs.get_conditions_metadata(conditions)
        conditions = shuffled_control_crfs.get_conditions_metadata(conditions)
        setup_logging(**conditions[sys.argv[1]])

        gridsearch_train_crfs.setup_exec_train_model(conditions)
        # Create bare-bones shuffle folder and create shuffled datasets
        shuffled_control_crfs.setup_shuffle_model(conditions)
        # Wait for train CRF to be done
        # Run merge and save_best, grabbing best params
        for cond in conditions.values():
            cond['to_test'] = gridsearch_train_crfs.test_train_CRFs
            cond['to_run'] = merge_and_get_parameters
        best_params = crf_util.wait_and_run(conditions)
        logger.info("Parameters for all conditions collected.\n")
        # create shuffle configs with best params (write and run write_configs_for_loopy.m)
        shuffled_control_crfs.create_shuffle_configs(conditions, best_params)
        # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
        for name, meta in conditions.items():
            meta['to_test'] = shuffled_control_crfs.simple_test_shuffle_datasets
            meta['to_run'] = shuffled_control_crfs.exec_shuffle_model
        crf_util.wait_and_run(conditions)
        # Wait for shuffle CRFs to be done and run merge and save_shuffle
        for meta in conditions.values():
            meta['to_test'] = shuffled_control_crfs.test_shuffle_CRFs
            meta['to_run'] = shuffled_control_crfs.exec_merge_shuffle_CRFs
        crf_util.wait_and_run(conditions)
        # Extract ensemble neuron IDs. Write to disk?
    else:
        raise TypeError("At least one condition name must be passed on the command line.")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
