#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Gridsearching on the Columbia habanero cluster.
"""
import time
import sys
import os
import stat
import subprocess
import logging

import crf_util
from gridsearch_train_crfs import GridsearchTrial


class GridsearchTrialHaba(GridsearchTrial):
    """Adaptation of GridsearchTrial to take advantage of Columbia's habanero cluster.

        Create an object and call .run() for basic usage.

    Attributes:
        email (str): Email for haba job submission. Set by settings file.
        email_jobs_threshold (int): Do not send email notifications for batch jobs larger
            than this, when email_notification is set to "num_jobs". Set by settings file.
        email_notification (str): Email notification parameter for haba scheduler. Set by
            settings file.
        group_id (str): Group ID for haba job submission. Set by settings file.
        username (str): Username for haba job submission. Set by settings file.
        haba_mem (str): Memory to request in MB. Set by HabaGridsearchOptions section of
            settings file. Passed verbatim to haba scheduler.
        haba_nodes (str): Passed verbatim to haba scheduler. Set by HabaGridsearchOptions
            section of settings file.
        haba_ppn (str): Passed verbatim to haba scheduler. Set by HabaGridsearchOptions
            section of settings file.
        haba_walltime (str): Walltime duration passed verbatim to habanero scheduler. Set by
            HabaGridsearchOptions section of settings file.
    """

    def __init__(
        self,
        condition_name,
        ini_fname="crf_parameters.ini",
        destination_path=None,
        logger=None,
    ):
        """
        Args:
            condition_name (str): String label for current trial.
            ini_fname (str, optional): Filepath to settings file. Defaults to
                "crf_parameters.ini" in current directory.
            destination_path (str, optional): Planned feature. No current use.
            logger (logging object, optional): Logger to use. Will produce its own by
                default.
        """
        if logger is None:
            logger = logging.getLogger("top." + __name__)
            logger.setLevel(logging.DEBUG)
        super().__init__(condition_name, ini_fname, destination_path, logger)

    def _parse_settings(self):
        """Override of Workflow. Calls super and collects any specificly required settings.
        """
        super()._parse_settings()
        parameters = crf_util.get_section_options(
            "HabaOptions", parser=self._parser, int_options=["email_jobs_threshold"]
        )
        self.username = parameters["username"]
        self.group_id = parameters["group_id"]
        self.email = parameters["email"]
        self.email_notification = parameters["email_notification"]
        self.email_jobs_threshold = parameters["email_jobs_threshold"]
        parameters.update(
            crf_util.get_section_options("HabaGridsearchOptions", parser=self._parser)
        )
        self.haba_nodes = parameters["haba_grid_nodes"]
        self.haba_ppn = parameters["haba_grid_ppn"]
        self.haba_walltime = parameters["haba_grid_walltime"]
        self.haba_mem = parameters["haba_grid_mem"]

    def test_train_CRFs(self):
        filebase = os.path.join(self._experiment, "results", "result")
        num_jobs = 1
        for param in [self.S_LAMBDAS, self.DENSITIES, self.P_LAMBDAS]:
            num_jobs *= param["num_points"] if param["parallize"] else 1
        return crf_util.get_max_job_done(filebase) >= num_jobs

    def _create_haba_config_sh(self):
        # TODO: Absolute paths
        # Expect to be in the experiment folder already when writing this
        num_jobs = 1
        for param in [self.S_LAMBDAS, self.DENSITIES, self.P_LAMBDAS]:
            num_jobs *= param["num_points"] if param["parallize"] else 1
        fname = "haba_config.sh"
        with open(fname, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#haba_config.sh\n\n")
            f.write("#Torque script to run Matlab program\n")
            f.write("\n#Torque directives\n")
            f.write("#SBATCH -J Gridsearch_{}\n".format(self._experiment))
            f.write("#SBATCH -A {}\n".format(self.group_id))
            f.write("#SBATCH -N {}\n".format(self.haba_nodes))
            f.write("#SBATCH --mem={}mb\n".format(self.haba_mem))
            f.write("#SBATCH --time={}\n".format(self.haba_walltime))
            if self.email_notification == "num_jobs":
                # Reduce email notifications for greater numbers of jobs
                if num_jobs == 1:
                    f.write("#SBATCH --mail-type=ALL\n")
                elif num_jobs <= self.email_jobs_threshold:
                    f.write("#SBATCH --mail-type=END\n")
                else:
                    f.write("#SBATCH --mail-type=FAIL\n")
            else:
                # Use email_notification setting verbatim
                f.write("#SBATCH --mail-type={}\n".format(self.email_notification))
            f.write("#SBATCH --mail-user={}\n".format(self.email))
            f.write("#SBATCH --array=1-{}\n".format(int(num_jobs)))
            # TODO: Use actual self.working_dir
            working_dir = os.path.join(self._expt_dir, self._experiment)
            f.write("#SBATCH --output={}/haba_logs/slurm-%A_%a.out\n".format(working_dir))
            f.write("#SBATCH --error={}/haba_logs/slurm-%A_%a.error\n".format(working_dir))
            f.write("cd {}\n".format(self.source_dir))
            f.write(
                "./run.sh {0} $SLURM_ARRAY_TASK_ID > expt/{0}/job_logs/matoutfile.$SLURM_ARRAY_TASK_ID\n".format(
                    self._experiment
                )
            )
            f.write("#End of script\n")
        f.closed
        # make sure file is executable:
        os.chmod(
            fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode
        )
        self._logger.info("Created " + fname + ".")

    def _create_start_jobs_sh(self, target="haba_config.sh"):
        # TODO: Absolute paths
        # Expect to be in the experiment folder already when writing this
        fname = "start_jobs.sh"
        with open(fname, "w") as f:
            # Clear out any basic remainders from previous runs
            f.write("rm -f ./results/result*.mat\n")
            f.write("rm -f ./haba_logs/*\n")
            f.write("rm -f ./job_logs/*\n")
            f.write(
                "cd ../.. && sbatch {}\n".format(
                    os.path.join("expt", self._experiment, target)
                )
            )
        f.closed
        # make sure file is executable
        os.chmod(
            fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode
        )
        self._logger.info("Created " + os.path.join(self._experiment, fname) + ".")

    def start_gridsearch_jobs(self):
        """Override of GridsearchTrial

        Raises:
            RuntimeError: If start_jobs.sh fails. Typically due to not being on habanero cluster.
        """
        self._create_haba_config_sh()
        self._create_start_jobs_sh()

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the habanero cluster? Job submission failed.")
            raise RuntimeError(
                "Received non-zero return code: {}".format(process_results)
            )
        self._logger.info("Training job(s) submitted.")

    def run(self):
        """Run a gridsearch trial.
        """
        # Update logging if module is invoked from the command line
        if __name__ == "__main__":
            # Create stdout log handler
            self._logger.addHandler(crf_util.get_StreamHandler(self.verbosity))
            self._logger.debug("Logging stream handler to sys.stdout added.")
        super().run()


if __name__ == "__main__":
    start_time = time.time()
    try:
        condition = sys.argv[1]
    except IndexError:
        raise TypeError("A condition name must be passed on the command line.")

    logger = logging.getLogger("top")
    logger.setLevel(logging.DEBUG)

    if len(sys.argv) > 2:
        gridsearch = GridsearchTrialHaba(condition, ini_fname=sys.argv[2], logger=logger)
    else:
        gridsearch = GridsearchTrialHaba(condition, logger=logger)
    gridsearch.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
