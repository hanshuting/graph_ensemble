#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Shuffled Controls trial on the Columbia yeti cluster.
"""
import time
import sys
import os
import stat
import subprocess
import logging

import crf_util
from shuffled_control_crfs import ShuffledControlsTrial


class ShuffledControlsTrialYeti(ShuffledControlsTrial):
    """docstring for ShuffledControlsTrialYeti"""

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
            crf_util.get_section_options(
                "YetiGenerateShuffledOptions", parser=self._parser
            )
        )
        self.yeti_gen_sh_nodes = parameters["yeti_gen_sh_nodes"]
        self.yeti_gen_sh_ppn = parameters["yeti_gen_sh_ppn"]
        self.yeti_gen_sh_walltime = parameters["yeti_gen_sh_walltime"]
        self.yeti_gen_sh_mem = parameters["yeti_gen_sh_mem"]
        parameters.update(
            crf_util.get_section_options(
                "YetiShuffledControlsOptions", parser=self._parser
            )
        )
        self.yeti_sh_ctrl_nodes = parameters["yeti_sh_ctrl_nodes"]
        self.yeti_sh_ctrl_ppn = parameters["yeti_sh_ctrl_ppn"]
        self.yeti_sh_ctrl_walltime = parameters["yeti_sh_ctrl_walltime"]
        self.yeti_sh_ctrl_mem = parameters["yeti_sh_ctrl_mem"]

    def create_shuffles(self):
        self._write_shuffling_yeti_script()
        self._write_shuffling_submit_script()
        self._logger.info(
            "done writing {} yeti scripts\n".format(self.shuffle_experiment)
        )

        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.shuffle_experiment)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))

        shell_command = ".{}shuffle_start_job.sh".format(os.sep)
        self._logger.debug("About to run {}".format(shell_command))
        process_results = subprocess.run(shell_command, shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError(
                "Received non-zero return code: {}".format(process_results)
            )
        self._logger.info("Shuffled dataset creation job submitted.")

        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def _write_shuffling_yeti_script(self):
        # write yeti script
        filepath = os.path.join(self.shuffle_experiment, "shuffle_yeti_config.sh")
        self._logger.debug("writing file: {}".format(filepath))
        with open(filepath, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#shuffle_yeti_config.sh\n")
            f.write(
                "#PBS -N Create_shuffled_dataset_{}\n".format(self.shuffle_experiment)
            )
            f.write("#PBS -W group_list={}\n".format(self.group_id))
            f.write(
                "#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
                    self.yeti_gen_sh_nodes,
                    self.yeti_gen_sh_ppn,
                    self.yeti_gen_sh_walltime,
                    self.yeti_gen_sh_mem,
                )
            )
            if self.email_notification == "num_jobs":
                # Always only 1 job, fully notify
                f.write("#PBS -m abe\n")
            else:
                # Use email_notification setting verbatim
                f.write("#PBS -m {}\n".format(self.email_notification))
            f.write("#PBS -M {}\n".format(self.email))
            f.write("#set output and error directories (SSCC example here)\n")
            log_folder = "{}/".format(
                os.path.join(self.shuffle_experiment_dir, "yeti_logs")
            )
            os.makedirs(os.path.expanduser(log_folder), exist_ok=True)
            f.write("#PBS -o localhost:{}\n".format(log_folder))
            f.write("#PBS -e localhost:{}\n".format(log_folder))
            f.write(
                "#Command below is to execute Matlab code for Job Array (Example 4) "
                + "so that each part writes own output\n"
            )
            f.write(
                'matlab -nodesktop -nodisplay -r "dbclear all;'
                + " addpath(genpath('{}'));".format(os.path.join(self.source_dir))
                + "gn_shuff_data('{}', '{}', {});".format(
                    os.path.join(self.data_dir, self.data_file),
                    self.shuffle_save_dir,
                    self.num_shuffle,
                )
                + ' exit"\n'
            )
            f.write("#End of script\n")
        f.closed
        # Set executable permissions
        os.chmod(
            filepath,
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode,
        )
        self._logger.info("done writing {}".format(filepath))

    def _write_shuffling_submit_script(self):
        # write submit script
        filepath = os.path.join(self.shuffle_experiment, "shuffle_start_job.sh")
        self._logger.debug("writing file: {}".format(filepath))
        with open(filepath, "w") as f:
            f.write("qsub shuffle_yeti_config.sh\n")
        f.closed
        # Set executable permissions
        os.chmod(
            filepath,
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode,
        )
        self._logger.info("done writing {}".format(filepath))

    def shuffle_training_prep(self):
        self._create_controls_yeti_config_sh()
        self._create_start_jobs_sh("controls_yeti_config.sh")

    def _create_controls_yeti_config_sh(self):
        # Expect to be in the experiment folder already when writing this
        # TODO: Absolute path from source_directory
        fname = "controls_yeti_config.sh"
        with open(fname, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#{}\n\n".format(fname))
            f.write("#Torque script to run Matlab program\n")
            f.write("\n#Torque directives\n")
            f.write("#PBS -N Shuffled_CRFs_{}\n".format(self.shuffle_experiment))
            f.write("#PBS -W group_list={}\n".format(self.group_id))
            f.write(
                "#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
                    self.yeti_sh_ctrl_nodes,
                    self.yeti_sh_ctrl_ppn,
                    self.yeti_sh_ctrl_walltime,
                    self.yeti_sh_ctrl_mem,
                )
            )
            num_jobs = int(self.num_shuffle)
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
            working_dir = os.path.expanduser(self.shuffle_experiment_dir)
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
                    self.shuffle_experiment
                )
            )
            f.write("#End of script\n")
        f.closed
        # make sure file is executable:
        os.chmod(
            fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode
        )
        self._logger.info("Created " + fname + ".")

    def _create_start_jobs_sh(self, target):
        # Expect to be in the experiment folder already when writing this
        fname = "start_jobs.sh"
        with open(fname, "w") as f:
            # Clear out any basic remainders from previous runs
            f.write("rm -f ./results/result*.mat\n")
            f.write("rm -f ./yeti_logs/*\n")
            f.write("rm -f ./job_logs/*\n")
            f.write(
                "cd ../.. && qsub {}\n".format(
                    os.path.join("expt", self.shuffle_experiment, target)
                )
            )
        f.closed
        # make sure file is executable
        os.chmod(
            fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode
        )
        self._logger.info("Created " + os.path.join(self.shuffle_experiment, fname) + ".")

    def train_controls(self):
        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self.shuffle_experiment_dir)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the yeti cluster? Job submission failed.")
            raise RuntimeError(
                "Received non-zero return code: {}".format(process_results)
            )
        self._logger.info("Shuffle jobs submitted.")
        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

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
        shuffled_trial = ShuffledControlsTrialYeti(
            condition, ini_fname=sys.argv[2], logger=logger
        )
    else:
        shuffled_trial = ShuffledControlsTrialYeti(condition, logger=logger)
    shuffled_trial.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
