import time
import sys
import os
import shlex
import subprocess
import configparser

import logging

logger = logging.getLogger("top." + __name__)
logger.setLevel(logging.DEBUG)


def loglevel_from_verbosity(verbosity):
    """Convert verbosity setting to logging level.

    Args:
        verbosity (int): Setting value. Meaningful values only between 0 and 5.
    """
    return max([logging.CRITICAL - (verbosity * 10), logging.DEBUG])


def get_FileHandler(log_fname, debug_filelogging=False, overwrite=True):
    """Setup logging to a log file on disk.

    Args:
        log_fname (str): File path at which to create the log file.
        debug_filelogging (bool, optional): If True, emit maximal messages to log. Default False,
            which emits one step down, omitting many large log records mostly used for debugging.
        overwrite (bool, optional): Determines whether to overwrite any preexisting log file.
            Default True.

    Returns:
        logging.FileHandler: The set up handler.
    """
    logfile_handler = logging.FileHandler(log_fname, mode="w" if overwrite else "a")
    logfile_handler.setLevel(logging.DEBUG if debug_filelogging else logging.INFO)
    logfile_format = logging.Formatter(
        "%(asctime)s - %(levelname)s@%(name)s: %(message)s"
    )
    logfile_handler.setFormatter(logfile_format)
    return logfile_handler


def get_StreamHandler(verbosity, stream=sys.stdout):
    """Setup logging to a stream, typically for viewing messages on screen.

    Args:
        verbosity (int): Sets how many messages to emit. Higher values produce more messages.
        stream (IO stream, optional): Where messages are emitted. Default to standard output.

    Returns:
        logging.StreamHandler: The set up handler.
    """
    stream_handler = logging.StreamHandler(stream=stream)
    stream_handler.setLevel(loglevel_from_verbosity(verbosity))
    return stream_handler


def get_raw_configparser(fname="crf_parameters.ini"):
    config = configparser.ConfigParser()
    logger.info("Loading settings file {}.".format(fname))
    config.read(fname)
    return config


def get_GridsearchOptions(parser=None, fname="crf_parameters.ini"):
    if parser is None:
        parser = get_raw_configparser(fname)
    GridsearchOptions = {
        section: {} for section in ["S_LAMBDAS", "DENSITIES", "P_LAMBDAS"]
    }
    for section, parameters in GridsearchOptions.items():
        parameters["parallize"] = parser.getboolean(section, "parallize")
        parameters["num_points"] = parser.getint(section, "num_points")
        parameters["min"] = parser.getfloat(section, "min")
        parameters["max"] = parser.getfloat(section, "max")
    return GridsearchOptions


def get_GeneralOptions(parser=None, fname="crf_parameters.ini"):
    if parser is None:
        parser = get_raw_configparser(fname)
    int_options = ["time_span", "verbosity"]
    bool_options_and_defaults = [
        ("debug_filelogging", False),
        ("no_same_neuron_edges", True),
    ]
    GeneralOptions = get_section_options(
        "GeneralOptions",
        parser=parser,
        int_options=int_options,
        bool_options_and_defaults=bool_options_and_defaults,
    )
    expanded_source = os.path.expanduser(GeneralOptions["source_directory"])
    if expanded_source != GeneralOptions["source_directory"]:
        logger.debug("Provided source_directory expanded to {}".format(expanded_source))
        GeneralOptions["source_directory"] = expanded_source
    expanded_data = os.path.expanduser(GeneralOptions["data_directory"])
    if expanded_data != GeneralOptions["data_directory"]:
        logger.debug("Provided data_directory expanded to {}".format(expanded_data))
        GeneralOptions["data_directory"] = expanded_data
    return GeneralOptions


def get_section_options(
    section,
    parser=None,
    fname="crf_parameters.ini",
    int_options=[],
    bool_options_and_defaults=[],
    float_options=[],
):
    """Pull settings into a dict, allowing for type specification.

    Args:
        section (str): Options file section name
        int_options (list of str, optional): Names of settings to coerce to int.
        bool_options_and_defaults (list of (str, bool), optional): Names of settings to coerce to
            bool and defaults to resort to if coersion fails. Ex. [("opt1", True), ("opt2", False)]
        float_options (list of str, optional): Names of settings to coerce to float.
        parser (ConfigParser, optional): Pre-existing parser to use.
        fname (str, optional): Filepath to options file to read if no parser is provided.

    Returns:
        dict: {option_name: value}
    """
    if parser is None:
        parser = get_raw_configparser(fname)
    section_options = {name: option for name, option in parser.items(section)}
    # Reread settings we expect to be non-string data types with specific getter
    for int_option in int_options:
        section_options[int_option] = parser.getint(section, int_option)
    for float_option in float_options:
        section_options[float_option] = parser.getfloat(section, float_option)
    for bool_option, option_default in bool_options_and_defaults:
        try:
            section_options[bool_option] = parser.getboolean(section, bool_option)
        except ValueError:
            section_options[bool_option] = option_default
    return section_options


def run_command(scommand, shell=False):
    logger.debug("About to run command:\n{}".format(scommand))
    if shell:
        process_results = subprocess.run(scommand, shell=True)
    else:
        sargs = shlex.split(scommand)
        process_results = subprocess.run(sargs)
    if process_results.returncode:
        raise RuntimeError("Received non-zero return code: {}".format(process_results))
    return process_results


def run_matlab_command(scommand, add_path=""):
    """Summary

    Args:
        scommand (str): Matlab command. Should end with a ";" or "," character.
        add_path (str, optional): If provided, this folder and all subfolders will be added to the
          Matlab path before running scommand.

    Returns:
        TYPE: Description
    """
    return run_command(
        'matlab -nodesktop -nodisplay -nosplash -r "'
        + ("addpath(genpath('{}')); ".format(add_path) if add_path else "")
        + scommand
        + 'exit"'
    )


def wait_and_run(to_test, to_run, wait_seconds=5):
    """Execute specified function after the corresponding test passes, pausing between tries.

    Args:
        to_test (() -> bool)): Parameter-less function that returns a bool value indicating
            when to execute to_run.
        to_run (() -> any)): Parameter-less function to run once to_test passes. Return value
            is passed as is.
        wait_seconds (float, optional): Number of seconds to wait per 'to_test' iterations.

    Returns:
        any: Return value from to_run.
    """
    logger.debug("Start waiting for: {}".format(to_test.__name__))
    num_waits = 0
    while not to_test():
        time.sleep(wait_seconds)
        num_waits += 1
        if (num_waits % 100) == 0:
            logger.info(
                "Waited {} sleep cycles so far testing {}".format(
                    num_waits, to_test.__name__
                )
            )
        elif (num_waits % 20) == 0:
            logger.debug(
                "Waited {} sleep cycles so far testing {}".format(
                    num_waits, to_test.__name__
                )
            )
    logger.info("to_test:{} passed.".format(to_test.__name__))
    logger.info("Now running {}.".format(to_run.__name__))
    return_val = to_run()
    logger.debug("return_val = {}".format(return_val))
    return return_val


def get_max_job_done(filebase, filesuffix=".mat"):
    filebase = os.path.expanduser(filebase)
    job = 1
    while os.path.exists("{}{}{}".format(filebase, job, filesuffix)):
        job += 1
    return job - 1
