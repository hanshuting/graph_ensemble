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
    return parameters


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


def create_start_jobs_sh(experiment):
    # Expect to be in the experiment folder already when writing this
    fname = "start_jobs.sh"
    with open(fname, 'w') as f:
        # Clear out any basic remainders from previous runs
        f.write("rm -f ./results/result*.mat\n")
        f.write("rm -f ./yeti_logs/*\n")
        f.write("rm -f ./job_logs/*\n")
        f.write("cd ../.. && qsub {}\n".format(os.path.join("expt", experiment, "yeti_config.sh")))
    f.closed

    # make sure file is executable
    os.chmod(fname, stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH | os.stat(fname).st_mode)
    logger.info("Created " + os.path.join(experiment, fname) + ".")


def start_gridsearch_jobs_yeti(params):
    create_yeti_config_sh(params)
    create_start_jobs_sh(params['experiment'])

    process_results = subprocess.run(".{}start_jobs.sh".format(os.sep), shell=True)
    return process_results
