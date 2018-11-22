"""Abstract superclass for model training workflows.
"""
import os
import logging

import crf_util


class Workflow(object):
    """Abstract superclass for all workflows.

    Attributes:
        condition_name (str): String label for current trial.
        data_dir (str): Path to dataset folder. Set by settings file.
        debug_filelogging (bool): Whether additional info will be emitted to log files.
            Set by settings file.
        destination_path (str): Planned feature. No current use.
        edges (str): Edge constraints setting. Set by settings file.
        experiment_group (str): Experiment group label for current trial. Set by settings
            file.
        ini_fname (str): Filepath to settings file.
        no_same_neuron_edges (bool): Ensure no edges between nodes for the same neuron,
            i.e. nodes across time_span for each neuron. Set by settings file.
        source_dir (str): Path to the root of this repo. Set in settings file.
        time_span (int): Maximum sequence length to consider. Set by settings file.
        verbosity (int): Determines degree of logging output to console. Set by settings
            file.
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
        self.condition_name = condition_name
        self.ini_fname = ini_fname
        if logger is None:
            logger = logging.getLogger("top." + __name__)
            logger.setLevel(logging.DEBUG)
        self._logger = logger
        if destination_path is None:
            self.destination_path = os.getcwd()
        else:
            self.destination_path = destination_path
        self._parse_settings()
        self._init_settings()

    def _parse_settings(self):
        """Reads in settings from settings file.
        """
        self._parser = crf_util.get_raw_configparser(fname=self.ini_fname)
        params = crf_util.get_GeneralOptions(parser=self._parser)
        self.experiment_group = params["experiment_name"]
        self.data_dir = params["data_directory"]
        self.source_dir = params["source_directory"]
        self.verbosity = params["verbosity"]
        self.debug_filelogging = params["debug_filelogging"]
        self.time_span = params["time_span"]
        self.edges = params["edges"]
        self.no_same_neuron_edges = params["no_same_neuron_edges"]

    def _init_settings(self):
        """Constructs internal handles.
        """
        # File name of dataset
        self._data_file = "{}_{}.mat".format(self.experiment_group, self.condition_name)
        self._experiment = "{}_{}".format(self.experiment_group, self.condition_name)
        # Path to expt folder
        self._expt_dir = os.path.join(self.source_dir, "expt")
