#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Class to gridsearch parameters for training a graphical model on neuronal
    spiketrain data.
"""
import time
import sys
import os
import logging

import crf_util
from workflow import Workflow


class GridsearchTrial(Workflow):
    """docstring for GridsearchTrial"""

    MODEL_TYPE = "loopy"
    # These parameters and their order must match best_parameters.txt.
    # See save_best_parameters.m for best_parameters.txt creation.
    PARAMS_TO_EXTRACT = ["s_lambda", "density", "p_lambda", "time_span"]

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
        """Override of Workflow
        """
        super()._parse_settings()
        params = crf_util.get_GridsearchOptions(parser=self._parser)
        self.S_LAMBDAS = params["S_LAMBDAS"]
        self.DENSITIES = params["DENSITIES"]
        self.P_LAMBDAS = params["P_LAMBDAS"]

    def _start_logfile(self):
        # TODO: Update to working_dir
        expt_dir = os.path.expanduser(self.expt_dir)
        log_fname = os.path.join(expt_dir, self.experiment, "gridsearch_train_crfs.log")
        logfile_handler = crf_util.get_FileHandler(
            log_fname, debug_filelogging=self.debug_filelogging
        )
        self._logger.addHandler(logfile_handler)
        self._logger.debug("Logging file handler to {} added.".format(log_fname))

    def _create_working_dir(self):
        # TODO: Use .destination_path
        self.working_dir = self.experiment
        self._logger.info("Creating working directory: {}".format(self.working_dir))
        os.makedirs(os.path.expanduser(self.working_dir))

    def _create_write_configs_for_loopy_m(self):
        # TODO: Fix fname path
        fname = os.path.join(self.experiment, "write_configs_for_loopy.m")
        # TODO: Just call create_configs directly
        with open(fname, "w") as f:
            f.write("create_configs( ...\n")
            f.write("    'datapath', '{}{}', ...\n".format(self.data_dir, self.data_file))
            f.write("    'experiment_name', '{}', ...\n".format(self.experiment))
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
            s_lambdas = self.S_LAMBDAS
            f.write(
                "    's_lambda_splits', {}, ...\n".format(
                    s_lambdas["num_points"] if s_lambdas["parallize"] else 1
                )
            )
            f.write(
                "    's_lambdas_per_split', {}, ...\n".format(
                    1 if s_lambdas["parallize"] else s_lambdas["num_points"]
                )
            )
            f.write("    's_lambda_min', {}, ...\n".format(s_lambdas["min"]))
            f.write("    's_lambda_max', {}, ...\n".format(s_lambdas["max"]))
            densities = self.DENSITIES
            f.write(
                "    'density_splits', {}, ...\n".format(
                    densities["num_points"] if densities["parallize"] else 1
                )
            )
            f.write(
                "    'densities_per_split', {}, ...\n".format(
                    1 if densities["parallize"] else densities["num_points"]
                )
            )
            f.write("    'density_min', {}, ...\n".format(densities["min"]))
            f.write("    'density_max', {}, ...\n".format(densities["max"]))
            p_lambdas = self.P_LAMBDAS
            f.write(
                "    'p_lambda_splits', {}, ...\n".format(
                    p_lambdas["num_points"] if p_lambdas["parallize"] else 1
                )
            )
            f.write(
                "    'p_lambdas_per_split', {}, ...\n".format(
                    1 if p_lambdas["parallize"] else p_lambdas["num_points"]
                )
            )
            f.write("    'p_lambda_min', {}, ...\n".format(p_lambdas["min"]))
            f.write("    'p_lambda_max', {}, ...\n".format(p_lambdas["max"]))
            f.write("    'edges', '{}', ...\n".format(self.edges.lower()))
            f.write(
                "    'no_same_neuron_edges', {}, ...\n".format(
                    str(self.no_same_neuron_edges).lower()
                )
            )
            f.write("    'time_span', {});\n".format(self.time_span))
        f.closed
        self._logger.info("done writing {}".format(fname))

    def start_gridsearch_jobs(self):
        if (
            self.S_LAMBDAS["parallize"]
            or self.DENSITIES["parallize"]
            or self.P_LAMBDAS["parallize"]
        ):
            self._logger.warning(
                "WARNING: ONLY config1.m WILL BE USED. "
                + "GridsearchTrial has no parallel execution capability, "
                + "but found parallize options set to True."
            )
        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.source_dir)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        shell_cmd = ".{}run.sh {} 1".format(os.sep, self.experiment)
        crf_util.run_command(shell_cmd, shell=True)
        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def setup_exec_train_model(self):
        """Mostly follows old create_script.pl.

            Expect to be in expt folder at start.

        Args:
            params (dict): Model parameters.
        """
        self._create_working_dir()
        self._start_logfile()
        self._create_write_configs_for_loopy_m()

        # Move into experiment folder
        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.experiment)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        crf_util.run_matlab_command(
            "write_configs_for_{},".format(self.MODEL_TYPE), add_path=self.source_dir
        )
        self._logger.info("\nTraining configs generated.")
        self.start_gridsearch_jobs()

        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def merge_save_train_models(self):
        results_path = os.path.join(self.experiment, "results")
        return crf_util.run_matlab_command(
            "save_best_params('{}'); ".format(results_path), add_path=self.source_dir
        )

    def get_best_parameters(self):
        """Get the best parameters from gridsearch results.

        Returns:
            dict: Best parameters from the gridsearch. Parameters stored as a dict with
                PARAMS_TO_EXTRACT as the keys.
        """
        best_params = {}
        # grab and return best params
        results_path = os.path.join(self.experiment, "results")
        with open(os.path.join(results_path, "best_parameters.txt"), "r") as f:
            for param in self.PARAMS_TO_EXTRACT:
                best_params[param] = float(f.readline())
        f.closed

        return best_params

    def test_train_CRFs(self):
        filebase = os.path.join(self.experiment, "results", "result")
        num_jobs = 1
        return crf_util.get_max_job_done(filebase) >= num_jobs

    def run(self):
        """Run a gridsearch trial.
        """
        if __name__ == "__main__":
            # Create stdout log handler if module is invoked from the command line
            self._logger.addHandler(crf_util.get_StreamHandler(self.verbosity))
            self._logger.debug("Logging stream handler to sys.stdout added.")

        self.setup_exec_train_model()
        # Wait for train CRF to be done and run merge and save_best
        crf_util.wait_and_run(self.test_train_CRFs, self.merge_save_train_models)
        best_params_path = os.path.join(
            self.expt_dir, self.experiment, "results", "best_parameters.txt"
        )
        self._logger.info(
            "Grid search complete. Best parameters in {}".format(best_params_path)
            + " in the following order:\n{}\n".format(self.PARAMS_TO_EXTRACT)
        )


if __name__ == "__main__":
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    logger = logging.getLogger("top")
    logger.setLevel(logging.DEBUG)

    if len(sys.argv) > 2:
        gridsearch = GridsearchTrial(condition, ini_fname=sys.argv[2], logger=logger)
    else:
        gridsearch = GridsearchTrial(condition, logger=logger)
    gridsearch.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
