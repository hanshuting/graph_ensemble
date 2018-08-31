import time
import sys
import os
import shlex
import subprocess
import configparser

import logging
logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
logger.setLevel(logging.INFO)


def loglevel_from_verbosity(verbosity):
    """Convert verbosity setting to logging level.

    Args:
        verbosity (int): Setting value. Meaningful values only between 0 and 5.
    """
    return max([logging.CRITICAL - (verbosity * 10), logging.DEBUG])


def get_raw_configparser(fname="crf_parameters.ini"):
    config = configparser.ConfigParser()
    config.read(fname)
    return config


def get_GridsearchOptions(parser=None, fname="crf_parameters.ini"):
    if parser is None:
        parser = get_raw_configparser(fname)
    GridsearchOptions = {section: {} for section in ["S_LAMBDAS", "DENSITIES", "P_LAMBDAS"]}
    for section, parameters in GridsearchOptions.items():
        parameters['parallize'] = parser.getboolean(section, 'parallize')
        parameters['num_points'] = parser.getint(section, 'num_points')
        parameters['min'] = parser.getfloat(section, 'min')
        parameters['max'] = parser.getfloat(section, 'max')
    return GridsearchOptions


def get_GeneralOptions(parser=None, fname="crf_parameters.ini"):
    if parser is None:
        parser = get_raw_configparser(fname)
    GeneralOptions = get_section_options("GeneralOptions", parser)
    # Reread settings we expect to be non-string data types with correct getter
    for int_option in ['time_span', 'num_shuffle', 'verbosity']:
        GeneralOptions[int_option] = parser.getint('GeneralOptions', int_option)
    for bool_option in ['debug_filelogging']:
        try:
            GeneralOptions[bool_option] = parser.getboolean('GeneralOptions', bool_option)
        except ValueError:
            GeneralOptions[bool_option] = False
    return GeneralOptions


def get_section_options(section, parser=None, fname="crf_parameters.ini"):
    if parser is None:
        parser = get_raw_configparser(fname)
    section_options = {name: option for name, option in parser.items(section)}
    return section_options


def get_user_parameters(fname="crf_parameters.ini"):
    config = get_raw_configparser(fname)
    parameters = {name: {} for name in config.sections()}
    for name, options in parameters.items():
        for option in config.options(name):
            parameters[name][option] = config.get(name, option)
    return parameters


def run_command(scommand):
    logger.debug("About to run:\n{}".format(scommand))
    sargs = shlex.split(scommand)
    process_results = subprocess.run(sargs)
    if process_results.returncode:
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    return process_results


def run_matlab_command(scommand, add_path=''):
    """Summary

    Args:
        scommand (str): Matlab command. Should end with a ";" or "," character.
        add_path (str, optional): If provided, this folder and all subfolders will be added to the
          Matlab path before running scommand.

    Returns:
        TYPE: Description
    """
    return run_command("matlab -nodesktop -nodisplay -nosplash -r \"" +
                       ("addpath(genpath('{}')); ".format(add_path) if add_path else '') +
                       scommand +
                       "exit\"")


def wait_and_run(conditions_to_check, wait_seconds=5):
    """Execute specified functions after their corresponding tests pass, pausing between tests.

    Args:
        conditions_to_check (dict of dicts): An item per waiting task and subsequent execution.
            Expects each top-level key to have an associated dict with at least:
                'to_test': a function that returns true when testing should conclude and execution
                    should begin.
                'to_run': the function to run once 'to_test' returns true. Ret
            The full dict of each top level key is passed as kwargs to its 'to_test' and 'to_run'.
        wait_seconds (float, optional): Number of seconds to wait per 'to_test' iterations.
    """
    return_vals = {}
    conditions_remaining = {name: None for name in conditions_to_check}
    logger.debug("Start waiting for\n{}".format(conditions_to_check))
    num_waits = 0
    while conditions_remaining:
        stop_checking = []
        for name in conditions_remaining:
            to_check = conditions_to_check[name]
            if to_check['to_test'](**to_check):
                logger.debug("{}['to_test'] passed.".format(name))
                # TODO: Parallize here so we can run but still continue to test others?
                return_vals[name] = to_check['to_run'](**to_check)
                logger.info("{} for {} completed.".format(to_check['to_run'].__name__, name))
                logger.debug("return_vals['{}'] = {}".format(name, return_vals[name]))
                stop_checking.append(name)
        for finished in stop_checking:
            del conditions_remaining[finished]

        time.sleep(wait_seconds)
        num_waits += 1
        if (num_waits % 100) == 0:
            logger.info("Waited for {} sleep cycles so far. Currently waiting for:\n{}".format(
                num_waits, conditions_to_check))
        elif (num_waits % 20) == 0:
            logger.debug("Waited for {} sleep cycles so far. Currently waiting for:\n{}".format(
                num_waits, conditions_to_check))

    logger.debug("Done waiting for {}.\n".format(conditions_to_check.keys()))
    # TODO: Returning all values together means the last test to pass blocks returing others.
    return return_vals


def get_max_job_done(filebase, filesuffix=".mat"):
    filebase = os.path.expanduser(filebase)
    job = 1
    while os.path.exists("{}{}{}".format(filebase, job, filesuffix)):
        job += 1
    return job - 1
