'''MIT License

Copyright (c) 2025 BioinfoMachineLearning

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.'''

# -------------------------------------------------------------------------------------------------------------------------------------
# Following code curated for PoseBench: (https://github.com/BioinfoMachineLearning/PSBench)
# -------------------------------------------------------------------------------------------------------------------------------------


import os, sys, argparse
import contextlib
import shutil
import tempfile
import time
from typing import Optional
from absl import logging


def die(msg):
    print(msg)
    sys.exit(1)


def check_dirs(params, keys, isdir=True):
    errmsg = ''
    for key in keys:
        dirpath = params[key]
        # print(f"{key}:{params[key]}")
        if isdir and not os.path.isdir(dirpath):
            errmsg = errmsg + '{}({})\n'.format(key, dirpath)

        if not isdir and not os.path.exists(dirpath):
            errmsg = errmsg + '{}({})\n'.format(key, dirpath)

    if len(errmsg) > 0:
        errmsg = 'Directories or files are not exist:\n' + errmsg
        raise argparse.ArgumentTypeError(errmsg)


def check_contents(params, keys):
    errmsg = ''
    for key in keys:
        name = params[key]
        if len(name) == 0:
            errmsg = errmsg + '{}\n'.format(key)

    if len(errmsg) > 0:
        errmsg = 'These contents are emply:\n' + errmsg
        raise argparse.ArgumentTypeError(msg)


def is_dir(dirname):
    """Checks if a path is an actual directory"""
    if not os.path.isdir(dirname):
        msg = "{0} is not a directory".format(dirname)
        raise argparse.ArgumentTypeError(msg)
    else:
        return dirname


def check_dir(dirname):
    return is_dir(dirname)


def is_file(filename):
    """Checks if a file is an invalid file"""
    if not os.path.exists(filename):
        msg = "{0} doesn't exist".format(filename)
        raise argparse.ArgumentTypeError(msg)
    else:
        return filename


def check_file(dirname):
    return is_file(dirname)


def makedir_if_not_exists(directory):
    os.makedirs(directory, exist_ok=True)

def run_command(cmd, log_file=None):
    flag = os.system(cmd)
    errmsg = 'Error occur while run: %s' % cmd
    if flag:
        print(errmsg)
        if log_file != None:
            write_str2file(log_file, errmsg)


def read_option_file(option_file):
    if not os.path.exists(option_file):
        die("Option file %s not exists." % option_file)
    params = {}
    for line in open(option_file):
        line = line.rstrip()
        if line.startswith('#'):
            continue
        tmp = line.split('=')
        if len(tmp) != 2:
            continue
        key = tmp[0].lstrip().rstrip()
        value = tmp[1].lstrip().rstrip()
        params[key] = value
    return params


def clean_dir(dir):
    if os.path.exists(dir):
        os.system(f'rm -rf {dir}')
    os.makedirs(dir)


def create_file(file):
    f = open(file, 'w')
    f.close()


def tmpdir_manager(base_dir):
    """Context manager that deletes a temporary directory on exit."""
    tmpdir = tempfile.mkdtemp(dir=base_dir)
    try:
        yield tmpdir
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def timing(msg: str):
    logging.info('Started %s', msg)
    tic = time.time()
    yield
    toc = time.time()
    logging.info('Finished %s in %.3f seconds', msg, toc - tic)

def extract_tmscore_file(tmscore_file):
    """
    Extracts TM-scores from a given US-align output file.

    Args:
        file_path (str): Path to the US-align output file.

    Returns:
        dict: A dictionary containing TM-scores normalized by Structure_1 and Structure_2.
              Example: {"Structure_1": 0.99536, "Structure_2": 0.99536}
    """
    
    try:
        with open(tmscore_file, 'r') as file:
            for line in file:
                # Check for TM-score normalized by Structure_1, ie: predicted structure
                if "TM-score=" in line and "normalized by length of Structure_1" in line:
                    tmscores_1 = float(line.split('=')[1].split()[0])
                
                # Check for TM-score normalized by Structure_2 ie: native structure
                elif "TM-score=" in line and "normalized by length of Structure_2" in line:
                    tmscores_2 = float(line.split('=')[1].split()[0])

        return(tmscores_2)
    
    except FileNotFoundError:
        print(f"File not found: {tmscore_file}")
        return "NA"
    except Exception as e:
        print(f"An error occurred: {e}")
        return "NA"