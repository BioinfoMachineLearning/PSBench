import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from util import is_file, is_dir, makedir_if_not_exists, clean_dir
import json

ema_exec = 'docker run --rm -v $(pwd):/home -v /bmlfast/pngkg/CASP16_DATASET_PAPER:/bmlfast/pngkg/CASP16_DATASET_PAPER/ -u $(id -u ${USER}):$(id -g ${USER}) registry.scicore.unibas.ch/schwede/openstructure:latest'

def run_command(inparams):
    mdl_path, trg_path, result_path = inparams

    cmd = [ema_exec, "compare-structures",
                       "-m", mdl_path, "-r", trg_path,
                       "-mf", "pdb", "-rf", "pdb", "-rna", "--out", result_path,
                       "--ics",
                       "--ips",
                       "--qs-score",
                       "--lddt",
                       "--local-lddt",
                       "--rigid-scores", 
                       "--patch-scores",
                       "--tm-score",
                       "--dockq",
                       "--ics-trimmed",
                       "--ips-trimmed"]

    cmd = ' '.join(cmd)
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True, help="Directory containing models for one target")
    parser.add_argument('--nativedir', type=is_file, required=True, help="Path to the native PDB file")
    parser.add_argument('--outdir', type=is_dir, required=True, help="Output directory")

    args = parser.parse_args()

    # Extract target name from the indir path
    targetname = os.path.basename(os.path.normpath(args.indir))

    # Define output directory for the target
    outdir = os.path.join(args.outdir, targetname)
    makedir_if_not_exists(outdir)

    # Native PDB file
    native_pdb = args.nativedir

    process_list = []

    for model in os.listdir(args.indir):
        infile = os.path.join(args.indir, model)
        outfile = os.path.join(outdir, f"{model}_filtered_out")

        # Check if the output file exists and contains valid data
        if os.path.exists(outfile):
            continue
            try:
                with open(outfile) as f:
                    data = f.read()
                    if 'qs_best' in data and 'tm_score' in data:
                        continue
            except Exception as e:
                print(f"Error reading {outfile}: {e}")

        process_list.append([args.indir + '/' + model, native_pdb, outfile])

    pool = Pool(processes=100)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
