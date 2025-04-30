from calendar import c
from curses import raw
import os, sys, argparse, time
from pydoc import doc
from multiprocessing import Pool
from tqdm import tqdm
import random
import numpy as np
from scipy.stats import pearsonr, spearmanr
import pandas as pd
from sklearn.metrics import mean_squared_error, roc_auc_score
from sklearn.preprocessing import MinMaxScaler

def process_target(target, group_id, args, native_df):
    """Processes a single target to compute correlation metrics."""
    
    prediction = os.path.join(args.indir, target)
    pred_df = pd.read_csv(prediction)
    pred_df['model'] = pred_df['model'] + ".pdb"

    common_models = set(pred_df['model'])
    common_native_df = native_df[native_df['model_name'].isin(common_models)]
    # print(common_native_df)

    scores_dict = {row['model_name']: float(row[args.native_score_field]) for _, row in common_native_df.iterrows()}
    true_tmscores = common_native_df[args.native_score_field].values
    
    if pred_df.empty:
        return ["0", "0", str(np.max(true_tmscores)), "0.5"]
    
    pred_df = pred_df.sort_values(by=[group_id], ascending=args.ascending).reset_index(drop=True)
    scores_filt, scores_true = [], []
    
    for i, row in pred_df.iterrows():
        model = row['model']
        if model in scores_dict:
            scores_filt.append(float(row[group_id]))
            scores_true.append(scores_dict[model])
    
    if not scores_filt:
        return ["0", "0", str(np.max(true_tmscores)),  "0.5"]
    
    corr = pearsonr(scores_filt, scores_true)[0]
    spear_corr = spearmanr(scores_filt, scores_true).statistic

    top1_model = pred_df.loc[0, 'model']
    if top1_model not in scores_dict:
        return ["0", "0", str(np.max(true_tmscores)), "0.5"]
    
    best_top1_tmscore = scores_dict[top1_model]
    loss = np.max(true_tmscores) - best_top1_tmscore
    roc_auc = max(0.5, roc_auc_score([int(x > np.quantile(true_tmscores, 0.75)) for x in scores_true], scores_filt))
    
    return [str(corr), str(spear_corr), str(loss), str(roc_auc)]

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--indir', type=str, required=True)
    parser.add_argument('--nativedir', type=str, required=True)
    parser.add_argument('--field', type=str, required=False)
    parser.add_argument('--ascending', type=bool, default=False, required=False)
    parser.add_argument('--native_score_field', type=str, default="usalign_score", required=False)
    args = parser.parse_args()
    
    scorefile = os.path.join(args.indir, os.listdir(args.indir)[0])
    group_ids = [args.field] if args.field else pd.read_csv(scorefile).columns[2:]
    print(group_ids)
    
    group_res = {}
    for group_id in group_ids:
        results = []
        for target in sorted(os.listdir(args.nativedir)):
            native_df = pd.read_csv(os.path.join(args.nativedir, target))
            targetname = target.replace('.csv', '')
            results.append(process_target(target, group_id, args, native_df))
        
        group_res[group_id] = {
            'corrs': [r[0] for r in results],
            'spear_corrs': [r[1] for r in results],
            'losses': [r[2] for r in results],
            'roc_aucs': [r[3] for r in results],
        }
    
    print('    '.join(group_res.keys()))
    targets = [target.rstrip('.csv') for target in sorted(os.listdir(args.nativedir))]
    print('\t'.join(targets))
    
    for i in range(len(os.listdir(args.nativedir))):
        print(' '.join(
            [val[i] for group_id in group_res for val in [group_res[group_id]['corrs'], group_res[group_id]['spear_corrs'],
                                                          group_res[group_id]['losses'], group_res[group_id]['roc_aucs']]]))

if __name__ == '__main__':
    main()

    
