import os, sys, argparse, time
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from util import is_file, is_dir, makedir_if_not_exists, clean_dir


def run_command(inparams):
    usalign_program, pdb1, pdb2, outfile = inparams
    cmd = f"{usalign_program} -ter 1 -tmscore 6 {pdb1} {pdb2} > {outfile} "
    print(cmd)
    os.system(cmd)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=is_dir, required=True, help="Directory containing models for one target")
    parser.add_argument('--nativedir', type=is_file, required=True, help="Path to the native PDB file")
    parser.add_argument('--outdir', type=is_dir, required=True, help="Output directory")
    parser.add_argument('--usalign_program', type=is_file, required=True, help="Path to the USAlign program")

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

        # Skip if the output file exists and has sufficient lines
        if os.path.exists(outfile) and len(open(outfile).readlines()) > 10:
            continue

        # Add the USAlign command to the process list
        process_list.append([args.usalign_program, infile, native_pdb, outfile])

    # Use multiprocessing to run USAlign commands
    pool = Pool(processes=180)
    results = pool.map(run_command, process_list)
    pool.close()
    pool.join()
