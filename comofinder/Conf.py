# -*- coding:utf-8 -*-
# project config
import os

from project.settings import UPLOAD_COMOFINDER_ROOT


# project data directory name
DATA_NAME = 'data'
# project output driectory name
OUTPUT_NAME = 'output'

RESULT_NAME = 'result'
# project log file driectory name
LOG_NAME = 'log'
# project version
VERSION = "1.0.1v"
# project usage help language: en, zh
LN = "en"
# project debug trigger
DEBUG = False
# motif class types: 0-miRNA, 1-TF, 2-gene
MOTIFTYPE = 3
# Maxrium Trials
MAXTRIALS = 100

# project current directory root path
PROJECT_DIR = UPLOAD_COMOFINDER_ROOT
# project data directory path
FILE_DIR = os.path.join(PROJECT_DIR, DATA_NAME)
# project running output directory path
OUTPUT_DIR = os.path.join(PROJECT_DIR, OUTPUT_NAME)

RESULT_DIR = os.path.join(PROJECT_DIR, RESULT_NAME)
# project logger file directory path
LOG_DIR = os.path.join(PROJECT_DIR, LOG_NAME)
