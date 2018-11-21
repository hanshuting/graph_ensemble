"""Abstract superclass for model training workflows.
"""
import os
import logging

import crf_util


class Workflow(object):
    """Abstract superclass for all workflows.

    Attributes:
        condition_name (str): Description
        ini_fname (str): Description
        logger (logger object): Description
    """

    def __init__(
        self,
        condition_name,
        ini_fname="crf_parameters.ini",
        destination_path=None,
        logger=None,
    ):
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
        self._expt_dir = os.path.join(self.source_dir, "expt")
