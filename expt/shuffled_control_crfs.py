#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import logging

import crf_util
from gridsearch_train_crfs import GridsearchTrial
import yeti_support
from workflow import Workflow


class ShuffledControlsTrial(Workflow):
    """docstring for ShuffledControlsTrial

    Attributes:
        condition_name (TYPE): Description
        ini_fname (TYPE): Description
        logger (TYPE): Description
        MODEL_TYPE (str): Description
    """

    MODEL_TYPE = "loopy"

    def __init__(
        self,
        condition_name,
        ini_fname="crf_parameters.ini",
        destination_path=None,
        logger=None,
    ):
        if logger is None:
            logger = logging.getLogger("top." + __name__)
            logger.setLevel(logging.DEBUG)
        super().__init__(condition_name, ini_fname, destination_path, logger)

    def _parse_settings(self):
        """Override of Workflow. Calls super and collects any specificly required settings.
        """
        super()._parse_settings()
        params = crf_util.get_GeneralOptions(parser=self._parser)
        self.cluster_architecture = params["cluster_architecture"]
        params = crf_util.get_section_options(
            "ShuffledControlsOptions", parser=self._parser, int_options=["num_shuffle"]
        )
        self.num_shuffle = params["num_shuffle"]

    def _init_settings(self):
        """Override of Workflow
        """
        super()._init_settings()
        self.shuffle_save_name = "shuffled_{}_{}".format(
            self.experiment_group, self.condition_name
        )
        self.shuffle_save_dir = os.path.join(
            os.path.expanduser(self.data_dir), "shuffled", self.experiment
        )
        self.shuffle_experiment = "shuffled_{}".format(self.experiment)
        self.shuffle_experiment_dir = os.path.join(
            self.source_dir, "expt", self.shuffle_experiment
        )
        # Update settings for cluster specified, if any
        if self.cluster_architecture == "yeti":
            self._logger.info("Yeti cluster architecture selected for shuffled controls.")
            yeti_support.get_yeti_shuff_metadata(self, parser=self._parser)

    def _start_logfile(self):
        # TODO: Update to working_dir
        log_fname = os.path.join(
            os.path.expanduser(self.shuffle_experiment_dir), "shuffled_control_crfs.log"
        )
        logfile_handler = crf_util.get_FileHandler(
            log_fname, debug_filelogging=self.debug_filelogging
        )
        self._logger.addHandler(logfile_handler)
        self._logger.debug("Logging file handler to {} added.".format(log_fname))

    def create_shuffles(self):
        matlab_cmd = "gn_shuff_data('{}', '{}', {});".format(
            os.path.join(self.data_dir, self.data_file),
            self.shuffle_save_dir,
            self.num_shuffle,
        )
        crf_util.run_matlab_command(matlab_cmd, add_path=self.source_dir)

    def setup_shuffle_model(self):
        self._logger.info(
            "Creating working directory: {}".format(self.shuffle_experiment)
        )
        os.makedirs(os.path.expanduser(self.shuffle_experiment))
        self._start_logfile()
        os.makedirs(os.path.expanduser(self.shuffle_save_dir), exist_ok=True)
        # TODO: Either clear pre-existing shuffled datasets, or skip regenerating any already there
        self.create_shuffles()

    def shuffle_training_prep(self):
        pass

    def create_shuffle_configs(self, best_params):
        fname = os.path.join(
            self.shuffle_experiment_dir, "write_shuffle_configs_for_loopy.m"
        )
        with open(fname, "w") as f:
            f.write("create_shuffle_configs( ...\n")
            f.write(
                "    'datapath', '{}.mat', ...\n".format(
                    os.path.join(self.shuffle_save_dir, self.shuffle_save_name)
                )
            )
            f.write("    'experiment_name', '{}', ...\n".format(self.shuffle_experiment))
            try:
                f.write("    'email_for_notifications', '{}', ...\n".format(self.email))
            except AttributeError:
                self._logger.debug("No notifications email setting provided. Skipping.")
            try:
                f.write("    'yeti_user', '{}', ...\n".format(self.username))
            except AttributeError:
                self._logger.debug("No yeti username provided. Skipping.")
            f.write("    'compute_true_logZ', false, ...\n")
            f.write("    'reweight_denominator', 'mean_degree', ...\n")
            # f.write("    's_lambda_splits', 1, ...\n")
            # f.write("    's_lambdas_per_split', 1, ...\n")
            f.write("    's_lambda_min', {}, ...\n".format(best_params["s_lambda"]))
            f.write("    's_lambda_max', {}, ...\n".format(best_params["s_lambda"]))
            # f.write("    'density_splits', 1, ...\n")
            # f.write("    'densities_per_split', 1, ...\n")
            f.write("    'density_min', {}, ...\n".format(best_params["density"]))
            f.write("    'density_max', {}, ...\n".format(best_params["density"]))
            # f.write("    'p_lambda_splits', 1, ...\n")
            # f.write("    'p_lambdas_per_split', 1, ...\n")
            f.write("    'p_lambda_min', {}, ...\n".format(best_params["p_lambda"]))
            f.write("    'p_lambda_max', {}, ...\n".format(best_params["p_lambda"]))
            f.write("    'edges', '{}', ...\n".format(self.edges.lower()))
            f.write(
                "    'no_same_neuron_edges', {}, ...\n".format(
                    str(self.no_same_neuron_edges).lower()
                )
            )
            f.write("    'time_span', {}, ...\n".format(best_params["time_span"]))
            f.write("    'num_shuffle', {});\n".format(self.num_shuffle))
        f.closed
        self._logger.info("done writing {}".format(fname))

        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.shuffle_experiment)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        crf_util.run_matlab_command(
            "write_shuffle_configs_for_{},".format(self.MODEL_TYPE),
            add_path=self.source_dir,
        )
        self.shuffle_training_prep()
        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def train_controls(self):
        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.source_dir)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        for cur_shuffle in range(1, self.num_shuffle + 1):
            scommand = ".{}run.sh {} {}".format(
                os.sep, self.shuffle_experiment, cur_shuffle
            )
            crf_util.run_command(scommand, shell=True)
            self._logger.info(
                "Trained shuffled control model {} out of {}.".format(
                    cur_shuffle, self.num_shuffle
                )
            )
        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def simple_test_shuffle_datasets(self):
        filebase = os.path.join(
            os.path.expanduser(self.shuffle_save_dir), self.shuffle_save_name + "_"
        )
        return crf_util.get_max_job_done(filebase) >= self.num_shuffle

    def test_shuffle_CRFs(self):
        filebase = os.path.join(self.shuffle_experiment, "results", "result")
        return crf_util.get_max_job_done(filebase) >= self.num_shuffle

    def exec_merge_shuffle_CRFs(self):
        results_path = os.path.join(self.shuffle_experiment, "results")
        crf_util.run_matlab_command(
            "save_and_merge_shuffled_models('{}'); ".format(results_path),
            add_path=self.source_dir,
        )
        self._logger.info("Shuffle models merged and saved.\n")

    def run(self):
        """Run a shuffled controls trial.
        """
        if __name__ == "__main__":
            # Create stdout log handler if module is invoked from the command line
            self._logger.addHandler(crf_util.get_StreamHandler(self.verbosity))
            self._logger.debug("Logging stream handler to sys.stdout added.")

        # Create bare-bones shuffle folder and create shuffled datasets
        self.setup_shuffle_model()
        # Get best params
        best_params = GridsearchTrial.get_best_parameters(
            os.path.join(self.expt_dir, self.experiment, "results")
        )
        self._logger.info("Parameters for {} collected.\n".format(self.condition_name))
        # create shuffle configs with best params (write and run write_configs_for_loopy.m)
        self.create_shuffle_configs(best_params)
        # Wait for all shuffled datasets to be created and run shuffle/start_jobs.sh
        crf_util.wait_and_run(self.simple_test_shuffle_datasets, self.train_controls)
        # Wait for shuffle CRFs to be done and run merge and save_shuffle
        crf_util.wait_and_run(self.test_shuffle_CRFs, self.exec_merge_shuffle_CRFs)
        # TODO: Extract ensemble neuron IDs and write to disk?


if __name__ == "__main__":
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    logger = logging.getLogger("top")
    logger.setLevel(logging.DEBUG)

    if len(sys.argv) > 2:
        shuffled_trial = ShuffledControlsTrial(
            condition, ini_fname=sys.argv[2], logger=logger
        )
    else:
        shuffled_trial = ShuffledControlsTrial(condition, logger=logger)
    shuffled_trial.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
