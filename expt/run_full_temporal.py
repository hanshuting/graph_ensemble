#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import crf_util
import logging
from gridsearch_train_crfs import GridsearchTrial
from shuffled_control_crfs import ShuffledControlsTrial


def main(gridsearch, shuffles, logger):
    gridsearch.setup_exec_train_model()
    # Create bare-bones shuffle folder and create shuffled datasets
    shuffles.setup_shuffle_model()
    # Wait for train CRF to be done
    # Run merge and save_best, grabbing best params
    crf_util.wait_and_run(gridsearch.test_train_CRFs, gridsearch.merge_save_train_models)
    best_params = GridsearchTrial.get_best_parameters(
        os.path.join(shuffles.expt_dir, shuffles.experiment, "results")
    )
    logger.info("Parameters for all conditions collected.\n")
    # create shuffle configs with best params (write and run write_configs_for_loopy.m)
    shuffles.create_shuffle_configs(best_params)
    # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
    crf_util.wait_and_run(shuffles.simple_test_shuffle_datasets, shuffles.train_controls)
    # Wait for shuffle CRFs to be done and run merge and save_shuffle
    crf_util.wait_and_run(shuffles.test_shuffle_CRFs, shuffles.exec_merge_shuffle_CRFs)
    # Extract ensemble neuron IDs. Write to disk?


if __name__ == "__main__":
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    logger = logging.getLogger("top")
    logger.setLevel(logging.DEBUG)
    init_log_handler = logging.StreamHandler(stream=sys.stdout)
    init_log_handler.setLevel(logging.INFO)
    logger.addHandler(init_log_handler)

    if len(sys.argv) > 2:
        ini_fname = sys.argv[2]
        logger.info("Using custom settings filename: {}".format(ini_fname))
        shuffles = ShuffledControlsTrial(condition, ini_fname=ini_fname)
        gridsearch = GridsearchTrial(condition, ini_fname=ini_fname)
    else:
        shuffles = ShuffledControlsTrial(condition)
        gridsearch = GridsearchTrial(condition)
    logger.addHandler(crf_util.get_StreamHandler(gridsearch.verbosity))
    logger.removeHandler(init_log_handler)
    main(gridsearch, shuffles, logger)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
