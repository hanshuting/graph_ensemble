"""Yeti specifics.
"""
import os
import stat
import subprocess
import crf_util

import logging
logger = logging.getLogger("top." + __name__)
logger.setLevel(logging.DEBUG)


def get_yeti_gs_metadata(fname="crf_parameters.ini"):
    parameters_parser = crf_util.get_raw_configparser(fname=fname)
    parameters = crf_util.get_section_options('YetiOptions', parser=parameters_parser)
    parameters.update(crf_util.get_section_options('YetiGridsearchOptions',
                                                   parser=parameters_parser))
    parameters['start_jobs'] = start_gridsearch_jobs_yeti
    parameters['test_gs_get_best_params'] = test_train_CRFs_yeti
    return parameters


def get_yeti_shuff_metadata(fname="crf_parameters.ini"):
    parameters_parser = crf_util.get_raw_configparser(fname=fname)
    parameters = crf_util.get_section_options('YetiOptions', parser=parameters_parser)
    parameters.update(crf_util.get_section_options('YetiGenerateShuffledOptions',
                                                   parser=parameters_parser))
    parameters.update(crf_util.get_section_options('YetiShuffledControlsOptions',
                                                   parser=parameters_parser))
    parameters["create_shuffles"] = create_shuffles_yeti
    parameters["shuffle_training_prep"] = shuff_controls_train_prep_yeti
    parameters["train_controls"] = exec_shuffle_model_yeti
    return parameters


def exec_shuffle_model_yeti(shuffle_experiment, **kwargs):
    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(shuffle_experiment)
    logger.debug("changed into dir: {}".format(os.getcwd()))

    process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
    if process_results.returncode:
        logger.critical("\nAre you on the yeti cluster? Job submission failed.")
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    logger.info("Shuffle jobs submitted.")

    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def shuff_controls_train_prep_yeti(params):
    create_controls_yeti_config_sh(params)
    create_start_jobs_sh(params['shuffle_experiment'], "controls_yeti_config.sh")


def create_controls_yeti_config_sh(params):
    # Expect to be in the experiment folder already when writing this
    # TODO: Absolute path from source_directory
    fname = "controls_yeti_config.sh"
    with open(fname, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#{}\n\n".format(fname))
        f.write("#Torque script to run Matlab program\n")

        f.write("\n#Torque directives\n")
        f.write("#PBS -N Shuffled_CRFs_{}\n".format(params['shuffle_experiment']))
        f.write("#PBS -W group_list={}\n".format(params['group_id']))
        f.write("#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
            params['yeti_sh_ctrl_nodes'], params['yeti_sh_ctrl_ppn'],
            params['yeti_sh_ctrl_walltime'], params['yeti_sh_ctrl_mem']))
        num_jobs = int(params['num_shuffle'])
        if params['email_notification'] == "num_jobs":
            # Reduce email notifications for greater numbers of jobs
            if num_jobs == 1:
                f.write("#PBS -m abe\n")
            elif num_jobs <= int(params['email_jobs_threshold']):
                f.write("#PBS -m ae\n")
            else:
                f.write("#PBS -m af\n")
        else:
            # Use email_notification setting verbatim
            f.write("#PBS -m {}\n".format(params['email_notification']))
        f.write("#PBS -M {}\n".format(params['email']))
        f.write("#PBS -t 1-{}\n".format(int(num_jobs)))

        working_dir = os.path.expanduser(params['shuffle_experiment_dir'])
        f.write("\n#set output and error directories (SSCC example here)\n")
        f.write("#PBS -o localhost:{}/yeti_logs/\n".format(working_dir))
        f.write("#PBS -e localhost:{}/yeti_logs/\n".format(working_dir))

        f.write("\n#Command below is to execute Matlab code for Job Array (Example 4) so that " +
                "each part writes own output\n")
        f.write("cd {}\n".format(params['source_directory']))
        f.write("./run.sh {0} $PBS_ARRAYID > expt/{0}/job_logs/matoutfile.$PBS_ARRAYID\n".format(
            params['shuffle_experiment']))
        f.write("#End of script\n")
    f.closed

    # make sure file is executable:
    os.chmod(fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode)
    logger.info("Created " + fname + ".")


def write_shuffling_yeti_script(params):
    # write yeti script
    filepath = os.path.join(params['shuffle_experiment'], "shuffle_yeti_config.sh")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#shuffle_yeti_config.sh\n")
        f.write("#PBS -N Create_shuffled_dataset_{}\n".format(params['shuffle_experiment']))
        f.write("#PBS -W group_list={}\n".format(params['group_id']))
        f.write("#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
            params['yeti_gen_sh_nodes'], params['yeti_gen_sh_ppn'],
            params['yeti_gen_sh_walltime'], params['yeti_gen_sh_mem']))
        if params['email_notification'] == "num_jobs":
            # Always only 1 job, fully notify
            f.write("#PBS -m abe\n")
        else:
            # Use email_notification setting verbatim
            f.write("#PBS -m {}\n".format(params['email_notification']))
        f.write("#PBS -M {}\n".format(params['email']))
        f.write("#set output and error directories (SSCC example here)\n")
        log_folder = "{}/".format(os.path.join(params['shuffle_experiment_dir'], "yeti_logs"))
        os.makedirs(os.path.expanduser(log_folder), exist_ok=True)
        f.write("#PBS -o localhost:{}\n".format(log_folder))
        f.write("#PBS -e localhost:{}\n".format(log_folder))
        f.write("#Command below is to execute Matlab code for Job Array (Example 4) " +
                "so that each part writes own output\n")
        f.write("matlab -nodesktop -nodisplay -r \"dbclear all;" +
                " addpath('{}');".format(os.path.join(params['shuffle_experiment_dir'])) +
                "gn_shuff_data('{}', '{}', {});".format(
                    os.path.join(params['data_directory'], params['data_file']),
                    params['shuffle_save_dir'],
                    params['num_shuffle']
                ) +
                " exit\"\n")
        f.write("#End of script\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}".format(filepath))


def write_shuffling_submit_script(experiment):
    # write submit script
    filepath = os.path.join(experiment, "shuffle_start_job.sh")
    logger.debug("writing file: {}".format(filepath))
    with open(filepath, 'w') as f:
        f.write("qsub shuffle_yeti_config.sh\n")
    f.closed
    # Set executable permissions
    os.chmod(filepath, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(filepath).st_mode)
    logger.info("done writing {}".format(filepath))


def create_shuffles_yeti(params):
    write_shuffling_yeti_script(params)
    write_shuffling_submit_script(params['shuffle_experiment'])
    logger.info("done writing {} yeti scripts\n".format(params['shuffle_experiment']))

    curr_dir = os.getcwd()
    logger.debug("curr_dir = {}.".format(curr_dir))
    os.chdir(params['shuffle_experiment'])
    logger.debug("changed into dir: {}".format(os.getcwd()))

    shell_command = ".{}shuffle_start_job.sh".format(os.sep)
    logger.debug("About to run {}".format(shell_command))
    process_results = subprocess.run(shell_command, shell=True)
    if process_results.returncode:
        logger.critical("\nAre you on the yeti cluster? Job submission failed.")
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    logger.info("Shuffled dataset creation job submitted.")

    os.chdir(curr_dir)
    logger.debug("changed back to dir: {}".format(os.getcwd()))


def test_train_CRFs_yeti(experiment, S_LAMBDAS, DENSITIES, P_LAMBDAS, **kwargs):
    filebase = os.path.join(experiment, "results", "result")
    num_jobs = 1
    for param in [S_LAMBDAS, DENSITIES, P_LAMBDAS]:
        num_jobs *= param['num_points'] if param['parallize'] else 1
    return crf_util.get_max_job_done(filebase) >= num_jobs


def create_yeti_config_sh(params):
    # Expect to be in the experiment folder already when writing this
    num_jobs = 1
    for param in [params['S_LAMBDAS'], params['DENSITIES'], params['P_LAMBDAS']]:
        num_jobs *= param['num_points'] if param['parallize'] else 1
    fname = "yeti_config.sh"
    with open(fname, 'w') as f:
        f.write("#!/bin/sh\n")
        f.write("#yeti_config.sh\n\n")
        f.write("#Torque script to run Matlab program\n")

        f.write("\n#Torque directives\n")
        f.write("#PBS -N Gridsearch_{}\n".format(params['experiment']))
        f.write("#PBS -W group_list={}\n".format(params['group_id']))
        f.write("#PBS -l nodes={}:ppn={},walltime={},mem={}mb\n".format(
            params['yeti_grid_nodes'], params['yeti_grid_ppn'],
            params['yeti_grid_walltime'], params['yeti_grid_mem']))
        if params['email_notification'] == "num_jobs":
            # Reduce email notifications for greater numbers of jobs
            if num_jobs == 1:
                f.write("#PBS -m abe\n")
            elif num_jobs <= int(params['email_jobs_threshold']):
                f.write("#PBS -m ae\n")
            else:
                f.write("#PBS -m af\n")
        else:
            # Use email_notification setting verbatim
            f.write("#PBS -m {}\n".format(params['email_notification']))
        f.write("#PBS -M {}\n".format(params['email']))
        f.write("#PBS -t 1-{}\n".format(int(num_jobs)))

        working_dir = os.path.join(params['expt_dir'], params['experiment'])
        f.write("\n#set output and error directories (SSCC example here)\n")
        f.write("#PBS -o localhost:{}/yeti_logs/\n".format(working_dir))
        f.write("#PBS -e localhost:{}/yeti_logs/\n".format(working_dir))

        f.write("\n#Command below is to execute Matlab code for Job Array (Example 4) so that " +
                "each part writes own output\n")
        f.write("cd {}\n".format(params['source_directory']))
        f.write("./run.sh {0} $PBS_ARRAYID > expt/{0}/job_logs/matoutfile.$PBS_ARRAYID\n".format(
            params['experiment']))
        f.write("#End of script\n")
    f.closed

    # make sure file is executable:
    os.chmod(fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode)
    logger.info("Created " + fname + ".")


def create_start_jobs_sh(experiment, target="yeti_config.sh"):
    # Expect to be in the experiment folder already when writing this
    fname = "start_jobs.sh"
    with open(fname, 'w') as f:
        # Clear out any basic remainders from previous runs
        f.write("rm -f ./results/result*.mat\n")
        f.write("rm -f ./yeti_logs/*\n")
        f.write("rm -f ./job_logs/*\n")
        f.write("cd ../.. && qsub {}\n".format(os.path.join("expt", experiment, target)))
    f.closed
    # make sure file is executable
    os.chmod(fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode)
    logger.info("Created " + os.path.join(experiment, fname) + ".")


def start_gridsearch_jobs_yeti(params):
    create_yeti_config_sh(params)
    create_start_jobs_sh(params['experiment'])

    process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
    if process_results.returncode:
        logger.critical("\nAre you on the yeti cluster? Job submission failed.")
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    logger.info("Training job(s) submitted.")
