#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Gridsearching on the Columbia yeti cluster.
"""
import time
import sys
import os
import stat
import subprocess
import logging

import crf_util
from gridsearch_train_crfs import GridsearchTrial


class GridsearchTrialYeti(GridsearchTrial):
    """docstring for GridsearchTrialYeti"""

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
        super()._parse_settings()
        parameters = crf_util.get_section_options("YetiOptions", parser=self._parser)
        self.username = parameters["username"]
        self.group_id = parameters["group_id"]
        self.email = parameters["email"]
        self.email_notification = parameters["email_notification"]
        self.email_jobs_threshold = parameters["email_jobs_threshold"]
        parameters.update(
            crf_util.get_section_options("YetiGridsearchOptions", parser=self._parser)
        )
        self.yeti_nodes = parameters["yeti_grid_nodes"]
        self.yeti_ppn = parameters["yeti_grid_ppn"]
        self.yeti_walltime = parameters["yeti_grid_walltime"]
        self.yeti_mem = parameters["yeti_grid_mem"]

    def test_train_CRFs(self):
        filebase = os.path.join(self._experiment, "results", "result")
        num_jobs = 1
        for param in [self.S_LAMBDAS, self.DENSITIES, self.P_LAMBDAS]:
            num_jobs *= param["num_points"] if param["parallize"] else 1
        return crf_util.get_max_job_done(filebase) >= num_jobs

    def _create_yeti_config_sh(self):
        # TODO: Absolute paths
        # Expect to be in the experiment folder already when writing this
        num_jobs = 1
        for param in [self.S_LAMBDAS, self.DENSITIES, self.P_LAMBDAS]:
            num_jobs *= param["num_points"] if param["parallize"] else 1
        fname = "yeti_config.sh"
        with open(fname, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#yeti_config.sh\n\n")
            f.write("#Torque script to run Matlab program\n")
            f.write("\n#Torque directives\n")
            f.write("#PBS -N Gridsearch_{}\n".format(self._experiment))
            f.write("#PBS -W group_list={}\n".format(self.group_id))
            f.write(
                "#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
                    self.yeti_nodes, self.yeti_ppn, self.yeti_walltime, self.yeti_mem
                )
            )
            if self.email_notification == "num_jobs":
                # Reduce email notifications for greater numbers of jobs
                if num_jobs == 1:
                    f.write("#PBS -m abe\n")
                elif num_jobs <= int(self.email_jobs_threshold):
                    f.write("#PBS -m ae\n")
                else:
                    f.write("#PBS -m af\n")
            else:
                # Use email_notification setting verbatim
                f.write("#PBS -m {}\n".format(self.email_notification))
            f.write("#PBS -M {}\n".format(self.email))
            f.write("#PBS -t 1-{}\n".format(int(num_jobs)))
            # TODO: Use actual self.working_dir
            working_dir = os.path.join(self._expt_dir, self._experiment)
            f.write("\n#set output and error directories (SSCC example here)\n")
            f.write("#PBS -o localhost:{}/yeti_logs/\n".format(working_dir))
            f.write("#PBS -e localhost:{}/yeti_logs/\n".format(working_dir))
            f.write(
                "\n#Command below is to execute Matlab code for Job Array (Example 4) so that "
                + "each part writes own output\n"
            )
            f.write("cd {}\n".format(self.source_dir))
            f.write(
                "./run.sh {0} $PBS_ARRAYID > expt/{0}/job_logs/matoutfile.$PBS_ARRAYID\n".format(
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

    def _create_start_jobs_sh(self, target="yeti_config.sh"):
        # TODO: Absolute paths
        # Expect to be in the experiment folder already when writing this
        fname = "start_jobs.sh"
        with open(fname, "w") as f:
            # Clear out any basic remainders from previous runs
            f.write("rm -f ./results/result*.mat\n")
            f.write("rm -f ./yeti_logs/*\n")
            f.write("rm -f ./job_logs/*\n")
            f.write(
                "cd ../.. && qsub {}\n".format(
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
            RuntimeError: If start_jobs.sh fails. Typically due to not being on yeti cluster.
        """
        self._create_yeti_config_sh()
        self._create_start_jobs_sh()

        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the yeti cluster? Job submission failed.")
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
        gridsearch = GridsearchTrialYeti(condition, ini_fname=sys.argv[2], logger=logger)
    else:
        gridsearch = GridsearchTrialYeti(condition, logger=logger)
    gridsearch.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
