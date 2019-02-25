#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""Shuffled Controls trial on the Columbia Haanero cluster.
"""
import time
import sys
import os
import stat
import subprocess
import logging

import crf_util
from shuffled_control_crfs import ShuffledControlsTrial


class ShuffledControlsTrialHaba(ShuffledControlsTrial):
    """Adaptation of ShuffledControlsTrial to take advantage of Columbia's Habanero cluster.

        Create an object and call .run() for basic usage.

    Attributes:
        email (str): Email for habanero job submission. Set by settings file.
        email_jobs_threshold (int): Do not send email notifications for batch jobs larger
            than this, when email_notification is set to "num_jobs". Set by settings file.
        email_notification (str): Email notification parameter for haba scheduler. Set by
            settings file.
        group_id (str): Group ID for haba job submission. Set by settings file.
        username (str): Username for haba job submission. Set by settings file.
        haba_gen_sh_mem (str): Memory to request in MB for creating shuffled datasets. Set
            by HabaGridsearchOptions section of settings file. Passed verbatim to haba
            scheduler.
        haba_gen_sh_nodes (str): For creating shuffled datasets. Passed verbatim to haba
            scheduler. Set by HabaGridsearchOptions section of settings file.
        haba_gen_sh_ppn (str): For creating shuffled datasets. Passed verbatim to haba
            scheduler. Set by HabaGridsearchOptions section of settings file.
        haba_gen_sh_walltime (str): Walltime duration for creating shuffled datasets
            passed verbatim to haba scheduler. Set by HabaGridsearchOptions section of
            settings file.
        haba_sh_ctrl_mem (str): Memory to request in MB for training models. Set by
            HabaGridsearchOptions section of settings file. Passed verbatim to haba
            scheduler.
        haba_sh_ctrl_nodes (str): For training models. Passed verbatim to haba scheduler.
            Set by HabaGridsearchOptions section of settings file.
        haba_sh_ctrl_ppn (str): For training models. Passed verbatim to haba scheduler.
            Set by HabaGridsearchOptions section of settings file.
        haba_sh_ctrl_walltime (str): Walltime duration for training models passed verbatim
            to haba scheduler. Set by HabaGridsearchOptions section of settings file.
    """

    def __init__(
        self,
        condition_name,
        ini_fname="crf_parameters.ini",
        destination_path=None,
        logger=None,
    ):
        """Summary

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
        parameters = crf_util.get_section_options("HabaOptions", parser=self._parser)
        self.username = parameters["username"]
        self.group_id = parameters["group_id"]
        self.email = parameters["email"]
        self.email_notification = parameters["email_notification"]
        self.email_jobs_threshold = parameters["email_jobs_threshold"]
        parameters.update(
            crf_util.get_section_options(
                "HabaGenerateShuffledOptions", parser=self._parser
            )
        )
        self.haba_gen_sh_nodes = parameters["haba_gen_sh_nodes"]
        self.haba_gen_sh_ppn = parameters["haba_gen_sh_ppn"]
        self.haba_gen_sh_walltime = parameters["haba_gen_sh_walltime"]
        self.haba_gen_sh_mem = parameters["haba_gen_sh_mem"]
        parameters.update(
            crf_util.get_section_options(
                "HabaShuffledControlsOptions", parser=self._parser
            )
        )
        self.haba_sh_ctrl_nodes = parameters["haba_sh_ctrl_nodes"]
        self.haba_sh_ctrl_ppn = parameters["haba_sh_ctrl_ppn"]
        self.haba_sh_ctrl_walltime = parameters["haba_sh_ctrl_walltime"]
        self.haba_sh_ctrl_mem = parameters["haba_sh_ctrl_mem"]

    def create_shuffles(self):
        self._write_shuffling_haba_script()
        self._write_shuffling_submit_script()
        self._logger.info(
            "done writing {} haba scripts\n".format(self._shuffle_experiment)
        )

        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self._shuffle_experiment)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))

        shell_command = ".{}shuffle_start_job.sh".format(os.sep)
        self._logger.debug("About to run {}".format(shell_command))
        process_results = subprocess.run(shell_command, shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the habanero cluster? Job submission failed.")
            raise RuntimeError(
                "Received non-zero return code: {}".format(process_results)
            )
        self._logger.info("Shuffled dataset creation job submitted.")

        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def _write_shuffling_haba_script(self):
        # write habanero script
        filepath = os.path.join(self._shuffle_experiment, "shuffle_haba_config.sh")
        self._logger.debug("writing file: {}".format(filepath))
        with open(filepath, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#shuffle_haba_config.sh\n")
            f.write(
                "#SBATCH --job-name=Create_shuffled_dataset_{}\n".format(self._shuffle_experiment)
            )
            f.write("#SBATCH --account={}".format(self.group_id))
            f.write("#SBATCH --nodes={}".format(self.haba_gen_sh_nodes))
            f.write("#SBATCH --time={}".format(self.haba_gen_sh_walltime))
            f.write("#SBATCH --mem {}".format(self.haba_gen_sh_mem))
            if self.email_notification == "num_jobs":
                # Always only 1 job, fully notify
                f.write("#SBATCH --mail-type=ALL\n")
            else:
                # Use email_notification setting verbatim
                f.write("#SBATCH --mail-type={}\n".format(self.email_notification))
            f.write("#SBATCH --mail-user={}\n".format(self.email))
            log_folder = "{}/".format(
                os.path.join(self._shuffle_experiment_dir, "haba_logs")
            )
            os.makedirs(os.path.expanduser(log_folder), exist_ok=True)
            f.write("#SBATCH --output={}slurm-%A_%a.out\n".format(log_folder))
            f.write("#SBATCH --error={}slurm-%A_%a.error\n".format(log_folder))
            f.write(
                'matlab -nodesktop -nodisplay -r "dbclear all;'
                + " addpath(genpath('{}'));".format(os.path.join(self.source_dir))
                + "gn_shuff_data('{}', '{}', {});".format(
                    os.path.join(self.data_dir, self._data_file),
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
        filepath = os.path.join(self._shuffle_experiment, "shuffle_start_job.sh")
        self._logger.debug("writing file: {}".format(filepath))
        with open(filepath, "w") as f:
            f.write("sbatch shuffle_haba_config.sh\n")
        f.closed
        # Set executable permissions
        os.chmod(
            filepath,
            stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode,
        )
        self._logger.info("done writing {}".format(filepath))

    def shuffle_training_prep(self):
        self._create_controls_haba_config_sh()
        self._create_start_jobs_sh("controls_haba_config.sh")

    def _create_controls_haba_config_sh(self):
        # Expect to be in the experiment folder already when writing this
        # TODO: Absolute path from source_directory
        fname = "controls_haba_config.sh"
        with open(fname, "w") as f:
            f.write("#!/bin/sh\n")
            f.write("#{}\n\n".format(fname))
            f.write("#SBATCH -J Shuffled_CRFs_{}\n".format(self._shuffle_experiment))
            f.write("#SBATCH -A {}\n".format(self.group_id))
            f.write("#SBATCH -t {}\n".format(self.haba_gen_sh_walltime))
            f.write("#SBATCH --mem={}\n".format(self.haba_gen_sh_mem))
            num_jobs = int(self.num_shuffle)
            if self.email_notification == "num_jobs":
                # Reduce email notifications for greater numbers of jobs
                if num_jobs == 1:
                    f.write("#SBATCH --mail-type=ALL\n")
                elif num_jobs <= int(self.email_jobs_threshold):
                    f.write("#SBATCH --mail-type=END\n")
                else:
                    f.write("#SBATCH --mail-type=FAIL\n")
            else:
                # Use email_notification setting verbatim
                f.write("#SBATCH --mail-type={}\n".format(self.email_notification))
            f.write("#SBATCH --mail-user={}\n".format(self.email))
            f.write("#SBATCH -a 1-{}\n".format(int(num_jobs)))
            working_dir = os.path.expanduser(self._shuffle_experiment_dir)
            f.write("#SBATCH --output={}/haba_logs/slurm-%A_%a.out\n".format(working_dir))
            f.write("#SBATCH --error={}/haba_logs/slurm-%A_%a.error\n".format(working_dir))
            f.write("cd {}\n".format(self.source_dir))
            f.write(
                "./run.sh {0} $SLURM_ARRAY_TASK_ID > expt/{0}/job_logs/matoutfile.$SLURM_ARRAY_TASK_ID\n".format(
                    self._shuffle_experiment
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
            f.write("rm -f ./haba_logs/*\n")
            f.write("rm -f ./job_logs/*\n")
            f.write(
                "cd ../.. && sbatch {}\n".format(
                    os.path.join("expt", self._shuffle_experiment, target)
                )
            )
        f.closed
        # make sure file is executable
        os.chmod(
            fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode
        )
        self._logger.info(
            "Created " + os.path.join(self._shuffle_experiment, fname) + "."
        )

    def train_controls(self):
        curr_dir = os.getcwd()
        self._logger.debug("curr_dir = {}.".format(curr_dir))
        os.chdir(self._shuffle_experiment_dir)
        self._logger.debug("changed into dir: {}".format(os.getcwd()))
        process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
        if process_results.returncode:
            self._logger.critical("\nAre you on the habanero cluster? Job submission failed.")
            raise RuntimeError(
                "Received non-zero return code: {}".format(process_results)
            )
        self._logger.info("Shuffle jobs submitted.")
        os.chdir(curr_dir)
        self._logger.debug("changed back to dir: {}".format(os.getcwd()))

    def run(self):
        """Run this shuffled control trial.
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
        shuffled_trial = ShuffledControlsTrialHaba(
            condition, ini_fname=sys.argv[2], logger=logger
        )
    else:
        shuffled_trial = ShuffledControlsTrialHaba(condition, logger=logger)
    shuffled_trial.run()

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
