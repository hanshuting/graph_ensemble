#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Expected to be run from the expt folder.
"""
import time
import sys
import os
import crf_util
import yeti_support

import logging

# *** start constants ***
MODEL_TYPE = "loopy"

# These parameters and their order must match best_parameters.txt.
# See save_best_parameters.m for best_parameters.txt creation.
PARAMS_TO_EXTRACT = ["s_lambda", "density", "p_lambda", "time_span"]
# *** end constants ***


def start_logfile(debug_filelogging, expt_dir, experiment, **_):
    expt_dir = os.path.expanduser(expt_dir)
    log_fname = os.path.join(expt_dir, experiment, "gridsearch_train_crfs.log")
    logfile_handler = crf_util.get_FileHandler(
        log_fname, debug_filelogging=debug_filelogging
    )
    logger.addHandler(logfile_handler)
    logger.debug("Logging file handler to {} added.".format(log_fname))


class GridsearchTrial(object):
    """docstring for GridsearchTrial"""

    def __init__(
        self,
        condition_name,
        ini_fname="crf_parameters.ini",
        destination_path=None,
        logger=None,
    ):
        # super(GridsearchTrial, self).__init__()
        self.condition_name = condition_name
        self.ini_fname = ini_fname
        if logger is None:
            logger = logging.getLogger("top." + __name__)
            logger.setLevel(logging.DEBUG)
        self.logger = logger

        self._init_settings()

    # def get_conditions_metadata(condition, ini_fname="crf_parameters.ini"):
    def _init_settings(self):
        """Reads in settings.

        Args:
            condition (str): Condition name.
            ini_fname (str, optional): Filepath of settings file to read.

        Returns:
            dict: Metadata for condition.
        """
        parameters_parser = crf_util.get_raw_configparser(fname=self.ini_fname)
        params = crf_util.get_GridsearchOptions(parser=parameters_parser)
        self.S_LAMBDAS = params["S_LAMBDAS"]
        self.DENSITIES = params["DENSITIES"]
        self.P_LAMBDAS = params["P_LAMBDAS"]
        params.update(crf_util.get_GeneralOptions(parser=parameters_parser))
        self.experiment_group = params["experiment_name"]
        self.data_dir = params["data_directory"]
        self.source_dir = params["source_directory"]
        self.verbosity = params["verbosity"]
        self.debug_filelogging = params["debug_filelogging"]
        self.cluster_architecture = params["cluster_architecture"]
        self.time_span = params["time_span"]
        self.edges = params["edges"]
        self.no_same_neuron_edges = params["no_same_neuron_edges"]
        self.data_file = "{}_{}.mat".format(params["experiment_name"], condition)
        self.experiment = "{}_{}_{}".format(
            params["experiment_name"], condition, MODEL_TYPE
        )
        self.expt_dir = os.path.join(params["source_directory"], "expt")

        self.start_jobs = start_gridsearch_jobs
        self.test_gs_get_best_params = test_train_CRFs
        # Update settings for cluster specified, if any
        if self.cluster_architecture == "yeti":
            self.logger.info("Yeti cluster architecture selected for gridsearch.")
            yeti_support.get_yeti_gs_metadata(self, fname=self.ini_fname)


def create_working_dir(params):
    logger.info("Creating working directory: {}".format(params["experiment"]))
    os.makedirs(os.path.expanduser(params["experiment"]))


def create_write_configs_for_loopy_m(params):
    fname = os.path.join(params["experiment"], "write_configs_for_loopy.m")
    # TODO: Just call create_configs directly
    with open(fname, "w") as f:
        f.write("create_configs( ...\n")
        f.write(
            "    'datapath', '{}{}', ...\n".format(
                params["data_directory"], params["data_file"]
            )
        )
        f.write("    'experiment_name', '{}', ...\n".format(params["experiment"]))
        try:
            f.write("    'email_for_notifications', '{}', ...\n".format(params["email"]))
        except KeyError:
            logger.debug("No notifications email setting provided. Skipping.")
        try:
            f.write("    'yeti_user', '{}', ...\n".format(params["username"]))
        except KeyError:
            logger.debug("No yeti username provided. Skipping.")
        f.write("    'compute_true_logZ', false, ...\n")
        f.write("    'reweight_denominator', 'mean_degree', ...\n")
        s_lambdas = params["S_LAMBDAS"]
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
        densities = params["DENSITIES"]
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
        p_lambdas = params["P_LAMBDAS"]
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
        f.write("    'edges', '{}', ...\n".format(params["edges"].lower()))
        f.write(
            "    'no_same_neuron_edges', {}, ...\n".format(
                str(params["no_same_neuron_edges"]).lower()
            )
        )
        f.write("    'time_span', {});\n".format(params["time_span"]))
    f.closed
    logger.info("done writing {}".format(fname))


def start_gridsearch_jobs(params):
    if (
        params["S_LAMBDAS"]["parallize"]
        or params["DENSITIES"]["parallize"]
        or params["P_LAMBDAS"]["parallize"]
    ):
        logger.warning(
            "WARNING: ONLY config1.m WILL BE USED. "
            + "Cluster architecture {} ".format(params["cluster_architecture"])
            + "has no parallel execution capability, "
            + "but found parallize options set to True."
        )
    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(params["source_directory"])
    logger.debug("changed into dir: {}".format(os.getcwd()))
    shell_cmd = ".{}run.sh {} 1".format(os.sep, params["experiment"])
    crf_util.run_command(shell_cmd, shell=True)
    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def setup_exec_train_model(params):
    """Mostly follows old create_script.pl.

        Expect to be in expt folder at start.

    Args:
        params (dict): Model parameters.
    """
    create_working_dir(params)
    start_logfile(**params)
    create_write_configs_for_loopy_m(params)

    # Move into experiment folder
    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(params["experiment"])
    logger.debug("changed into dir: {}".format(os.getcwd()))
    crf_util.run_matlab_command(
        "write_configs_for_{},".format(MODEL_TYPE), add_path=params["source_directory"]
    )
    logger.info("\nTraining configs generated.")

    # start_gridsearch_jobs by default
    params["start_jobs"](params)

    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def merge_save_train_models(experiment, source_directory, **kwargs):
    results_path = os.path.join(experiment, "results")
    return crf_util.run_matlab_command(
        "save_best_params('{}'); ".format(results_path), add_path=source_directory
    )


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
    with open(os.path.join(results_path, "best_parameters.txt"), "r") as f:
        for param in PARAMS_TO_EXTRACT:
            best_params[param] = float(f.readline())
    f.closed

    return best_params


def test_train_CRFs(experiment, **kwargs):
    filebase = os.path.join(experiment, "results", "result")
    num_jobs = 1
    return crf_util.get_max_job_done(filebase) >= num_jobs


def main(condition, ini_fname="crf_parameters.ini"):
    """Summary

    Args:
        condition (str): Condition name.
        ini_fname (str, optional): Filepath of settings file to read.
    """
    params = get_conditions_metadata(condition, ini_fname)

    # Update logging if module is invoked from the command line
    if __name__ == "__main__":
        # Assume top log position
        logger = logging.getLogger("top")
        logger.setLevel(logging.DEBUG)
        # Create stdout log handler
        verbosity = params["verbosity"]
        logger.addHandler(crf_util.get_StreamHandler(verbosity))
        logger.debug("Logging stream handler to sys.stdout added.")

    setup_exec_train_model(params)
    # Wait for train CRF to be done
    # Run merge and save_best
    params["to_test"] = params["test_gs_get_best_params"]
    params["to_run"] = merge_save_train_models
    crf_util.wait_and_run(params)
    best_params_path = os.path.join(
        params["expt_dir"], params["experiment"], "results", "best_parameters.txt"
    )
    logger.info(
        "Grid search complete. Best parameters in {}".format(best_params_path)
        + " in the following order:\n{}\n".format(PARAMS_TO_EXTRACT)
    )


if __name__ == "__main__":
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    if len(sys.argv) > 2:
        main(condition, sys.argv[2])
    else:
        main(condition)

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
