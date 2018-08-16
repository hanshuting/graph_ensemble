#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import crf_util
import gridsearch_train_crfs
import shuffled_control_crfs

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))


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
        verbosity = int(conditions[sys.argv[1]]['verbosity'])
        logger.setLevel(crf_util.loglevel_from_verbosity(verbosity))
        conditions = shuffled_control_crfs.get_conditions_metadata(conditions)
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
