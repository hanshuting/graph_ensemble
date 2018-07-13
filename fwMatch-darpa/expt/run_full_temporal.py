#!/usr/bin/python
# -*- coding: utf-8 -*-
import time
import sys
import os
import shutil
import subprocess

# *** START USER EDITABLE VARIABLES ***
EXPT_NAME = "temporal"
DATA_DIR = "~/data/"
USER = "jds2270"
EMAIL = "jds2270@columbia.edu"
# *** END USER EDITABLE VARIABLES ***


MODEL_TYPE = "loopy"
TRAIN_TEMPLATE_FOLDER_NAME = "{}_template".format(EXPT_NAME)
#SHUFFLE_TEMPLATE_FOLDER_NAME = # TODO


def check_templates():
    """Make template folders from master templates if they do not exist already
    """
    # Check train
    if not os.path.exists(TRAIN_TEMPLATE_FOLDER_NAME):
        shutil.copytree("temporal_template", TRAIN_TEMPLATE_FOLDER_NAME)

    # TODO: Check stuffle


def setup_train_working_dir(condition_names):
    for condition in condition_names:
        data_file = "{}_{}.mat".format(EXPT_NAME, condition)
        experiment = "{}_{}_{}".format(EXPT_NAME, condition, MODEL_TYPE)
        print("Copying {} to {}".format(TRAIN_TEMPLATE_FOLDER_NAME, experiment))
        shutil.copytree(TRAIN_TEMPLATE_FOLDER_NAME, experiment)

        with open("{}{}write_configs_for_loopy.m".format(experiment, os.sep), 'w') as f:
            f.write("create_config_files( ...\n")
            f.write("    'datapath', '{}{}', ...\n".format(DATA_DIR, data_file))
            f.write("    'experiment_name', '{}', ...\n".format(experiment))
            f.write("    'email_for_notifications', '{}', ...\n".format(EMAIL))
            f.write("    'yeti_user', '{}', ...\n".format(USER))
            f.write("    'compute_true_logZ', false, ...\n")
            f.write("    'reweight_denominator', 'mean_degree', ...\n")
            f.write("    's_lambda_splits', 6, ...\n")
            f.write("    's_lambdas_per_split', 1, ...\n")
            f.write("    's_lambda_min', 2e-03, ...\n")
            f.write("    's_lambda_max', 5e-01, ...\n")
            f.write("    'density_splits', 1, ...\n")
            f.write("    'densities_per_split', 6, ...\n")
            f.write("    'density_min', 0.1, ...\n")
            f.write("    'density_max', 0.3, ...\n")
            f.write("    'p_lambda_splits', 5, ...\n")
            f.write("    'p_lambdas_per_split', 1, ...\n")
            f.write("    'p_lambda_min', 1e+01, ...\n")
            f.write("    'p_lambda_max', 1e+04, ...\n")
            f.write("    'time_span', 2);\n")
        f.closed
        print("done writing write_configs_for_loopy.m\n")

        curr_dir = os.curdir()

if __name__ == '__main__':
    start_time = time.time()
    conditions = sys.argv[1:]
    if conditions:
        check_templates()
        setup_train_working_dir(conditions)
    else:
        raise TypeError("At least one condition name must be passed on the command line.\n")

    print("Total run time: {0:.2f} seconds".format(time.time() - start_time))
