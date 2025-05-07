# -------------------------------------------------------------------------------------------------------------------------------------
# Following code curated for PoseBench: (https://github.com/BioinfoMachineLearning/PSBench)
# -------------------------------------------------------------------------------------------------------------------------------------


import os
import pickle
import numpy as np
from tqdm import tqdm
import argparse
from util import is_dir,is_file,makedir_if_not_exists
from multiprocessing import Pool, cpu_count

def extract_scores(pdb_path, pkl_path, fasta_path, pae_cutoff=5.0):
    def parse_atm_record(line):
        return {
            'atm_name': line[12:16].strip(),
            'res_name': line[17:20].strip(),
            'chain': line[21],
            'x': float(line[30:38]),
            'y': float(line[38:46]),
            'z': float(line[46:54]),
        }

    def read_pdb(pdbfile):
        coords, cb_inds, ca_inds = {}, {}, {}
        with open(pdbfile) as file:
            for line in file:
                if 'ATOM' in line:
                    r = parse_atm_record(line)
                    coords.setdefault(r['chain'], []).append([r['x'], r['y'], r['z']])
                    ca_inds.setdefault(r['chain'], [])
                    cb_inds.setdefault(r['chain'], [])
                    idx = len(coords[r['chain']]) - 1
                    if r['atm_name'] == 'CA':
                        ca_inds[r['chain']].append(idx)
                    if r['atm_name'] == 'CB' or (r['atm_name'] == 'CA' and r['res_name'] == 'GLY'):
                        cb_inds[r['chain']].append(idx)
        return coords, cb_inds, ca_inds

    def read_plddt(plddt, ca_inds):
        plddt_chain = {}
        offset = 0
        for chain, indices in ca_inds.items():
            length = len(indices)
            plddt_chain[chain] = plddt[offset:offset+length]
            offset += length
        return plddt_chain

    def score_complex(coords, cb_inds, plddt_chain):
        chains = list(coords.keys())
        total_score = 0
        for i, c1 in enumerate(chains):
            cb1 = np.array(coords[c1])[cb_inds[c1]]
            plddt1 = plddt_chain[c1]
            l1 = len(cb1)
            for j, c2 in enumerate(chains):
                if j == i: continue
                cb2 = np.array(coords[c2])[cb_inds[c2]]
                plddt2 = plddt_chain[c2]
                mat = np.concatenate((cb1, cb2), axis=0)
                d = np.sqrt(np.sum((mat[:, None] - mat[None]) ** 2, axis=2))
                contact = d[:l1, l1:] <= 8
                if np.any(contact):
                    ids = np.argwhere(contact)
                    avg_plddt = np.mean(np.concatenate([plddt1[ids[:,0]], plddt2[ids[:,1]]]))
                    total_score += np.log10(len(ids) + 1) * avg_plddt
        return total_score, len(chains)

    def mpDockQ(score):
        L, x0, k, b = 0.827, 261.398, 0.036, 0.221
        return L / (1 + np.exp(-k*(score - x0))) + b

    def pDockQ(coords, cb_inds, plddt_chain, t=8):
        ch1, ch2 = list(cb_inds.keys())
        cb1, cb2 = np.array(coords[ch1])[cb_inds[ch1]], np.array(coords[ch2])[cb_inds[ch2]]
        d = np.sqrt(np.sum((cb1[:, None] - cb2[None])**2, axis=-1))
        ids = np.argwhere(d <= t)
        if len(ids) < 1: return 0.0
        avg_plddt = np.mean(np.concatenate([plddt_chain[ch1][np.unique(ids[:,0])], plddt_chain[ch2][np.unique(ids[:,1])]]))
        x = avg_plddt * np.log10(len(ids))
        return 0.724 / (1 + np.exp(-0.052 * (x - 152.611))) + 0.018

    def count_inter_pae(pae, seqs, cutoff):
        lens = [len(s) for s in seqs]
        offset = 0
        for l in lens:
            pae[offset:offset+l, offset:offset+l] = 50
            offset += l
        return np.sum(pae < cutoff)

    # Load data
    seqs = [line.rstrip('\n') for line in open(fasta_path) if line[0] != '>']
    # print(pkl_path)
    result = pickle.load(open(str(pkl_path), 'rb'))
    # print(pkl_path)
    iptm_ptm = float(result['ranking_confidence'])
    iptm = float(result['iptm'])
    pae = result['predicted_aligned_error']
    plDDT = result['plddt']

    coords, cb_inds, ca_inds = read_pdb(pdb_path)
    plddt_chain = read_plddt(plDDT, ca_inds)
    score, nchains = score_complex(coords, cb_inds, plddt_chain)
    mpdockq = mpDockQ(score) if nchains > 2 else pDockQ(coords, cb_inds, plddt_chain)
    num_inter_pae = count_inter_pae(pae, seqs, pae_cutoff)

    return {
        "iptm_ptm": iptm_ptm,
        "iptm": iptm,
        "num_inter_pae": num_inter_pae,
        "mpDockQ/pDockQ": mpdockq
    }

import pandas as pd
import os


def process_pair(args):
    """Process a single model by matching .pdb and .pkl files."""
    pdb_path, pkl_path, fasta_path = args
    model_name = os.path.basename(pdb_path)
    try:
        scores = extract_scores(pdb_path, pkl_path, fasta_path)
        return {'model_name': model_name, **scores}
    except Exception as e:
        print(f"Failed processing {model_name}: {e}")
        return None

def generate_af_features(fasta_path,pdb_dir, pkl_dir, output_csv, n_workers=None):
    os.makedirs(os.path.dirname(output_csv), exist_ok=True)
    n_workers = n_workers or cpu_count()

    pdb_files = [f for f in os.listdir(pdb_dir) if f.endswith('.pdb')]
    tasks = []

    for pdb_file in pdb_files:
        model_name = os.path.splitext(pdb_file)[0]
        pdb_path = os.path.join(pdb_dir, pdb_file)
        pkl_path = os.path.join(pkl_dir, model_name + '.pkl')

        if os.path.exists(pkl_path):
            tasks.append((pdb_path, pkl_path, fasta_path))
        else:
            print(f"Missing .pkl for {model_name}, skipping.")

    with Pool(n_workers) as pool:
        results = list(tqdm(pool.imap_unordered(process_pair, tasks), total=len(tasks)))

    results = [r for r in results if r is not None]
    df = pd.DataFrame(results)
    df.to_csv(output_csv, index=False)
    print(f"Saved results to {output_csv}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=is_file, required=True, help="Directory of Fasta file")
    parser.add_argument('--predicted_dir', type=is_dir, required=True, help="Directory of original model PDBs")
    parser.add_argument('--pkl_dir', type=is_dir, required=True, help="Directory of original model pkls")
    parser.add_argument('--outdir', required=True, help="Path to output csv")

    args = parser.parse_args()
    makedir_if_not_exists(args.outdir)
    targetname = os.path.basename(os.path.normpath(args.predicted_dir))
    outcsv_path = os.path.join(args.outdir,f"{targetname}_af_features.csv")
    generate_af_features(args.fasta,args.predicted_dir,args.pkl_dir,outcsv_path)
