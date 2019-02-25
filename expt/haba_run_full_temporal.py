#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Unified gridsearch and shuffled controls. Leverages asynchronous scheduling of Columbia
    Habanero cluster to optimally interleave subtasks of each trial type.

    Run from the command line with a condition name argument. Optionally, pass the name of
        the settings file to use as a second argument.
"""
import time
import sys
import crf_util
import logging
from haba_gridsearch_train_crfs import GridsearchTrialHaba
from haba_shuffled_control_crfs import ShuffledControlsTrialHaba
import run_full_temporal


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
        shuffles = ShuffledControlsTrialHaba(condition, ini_fname=ini_fname)
        gridsearch = GridsearchTrialHaba(condition, ini_fname=ini_fname)
    else:
        shuffles = ShuffledControlsTrialHaba(condition)
        gridsearch = GridsearchTrialHaba(condition)
    logger.addHandler(crf_util.get_StreamHandler(gridsearch.verbosity))
    logger.removeHandler(init_log_handler)
    run_full_temporal.main(gridsearch, shuffles, logger)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
