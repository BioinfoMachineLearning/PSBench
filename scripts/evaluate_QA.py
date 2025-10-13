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
# Following code curated for PSBench: (https://github.com/BioinfoMachineLearning/PSBench)
# -------------------------------------------------------------------------------------------------------------------------------------


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

def process_target(target, ema_method, args, native_df):
    """Processes a single target to compute correlation metrics."""
    
    prediction = os.path.join(args.input_dir, target + '.csv')
    pred_df = pd.read_csv(prediction)
    pred_df['model'] = pred_df['model'].apply(lambda x: x if x.endswith('.pdb') else f"{x}.pdb")
    # print(pred_df)
    common_models = set(pred_df['model'])
    common_native_df = native_df[native_df['model_name'].isin(common_models)]
    # print(common_native_df)

    scores_dict = {row['model_name']: float(row[args.true_score_field]) for _, row in common_native_df.iterrows()}
    true_tmscores = common_native_df[args.true_score_field].values
    
    if pred_df.empty:
        return ["0", "0", str(np.max(true_tmscores)), "0.5"]
    
    if ema_method not in pred_df.columns:
        return ["0", "0", str(np.max(true_tmscores)),  "0.5"]
        
    pred_df = pred_df.sort_values(by=[ema_method], ascending=False).reset_index(drop=True)
    scores_filt, scores_true = [], []
    
    for i, row in pred_df.iterrows():
        model = row['model']
        if model in scores_dict:
            scores_filt.append(float(row[ema_method]))
            scores_true.append(scores_dict[model])
    
    if not scores_filt:
        return ["0", "0", str(np.max(true_tmscores)),  "0.5"]
    
    if np.any(np.isnan(scores_filt)) or np.any(np.isnan(scores_true)):
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
    parser.add_argument('--input_dir', type=str, required=True)
    parser.add_argument('--native_dir', type=str, required=True)
    parser.add_argument('--ema_method', type=str, required=False)
    parser.add_argument('--true_score_field', type=str, default="tmscore_usalign", required=False)
    parser.add_argument('--outfile', type=str, default="evaluation_results.csv", help="Path to output CSV file")
    args = parser.parse_args()
    
    scorefile = os.path.join(args.input_dir, os.listdir(args.input_dir)[0])
    ema_methods = [args.ema_method] if args.ema_method else [column for column in pd.read_csv(scorefile).columns if column != 'model']
    print(f"EMA methods for evaluation: {', '.join(ema_methods)}")
    
    targets = [target.replace('.csv', '') for target in sorted(os.listdir(args.input_dir))]
    all_rows = []

    for i, target in enumerate(targets):
        native_df = pd.read_csv(os.path.join(args.native_dir, target + '_quality_scores.csv'))
        row = {"target": target}
        for ema_method in ema_methods:
            result = process_target(target, ema_method, args, native_df)
            row[f"{ema_method}_pearson"] = result[0]
            row[f"{ema_method}_spearman"] = result[1]
            row[f"{ema_method}_loss"] = result[2]
            row[f"{ema_method}_auroc"] = result[3]
        all_rows.append(row)

    # Create and save the DataFrame
    result_df = pd.DataFrame(all_rows)
    result_df.to_csv(args.outfile, index=False)
    print(f"\nResults written to {args.outfile}")

if __name__ == '__main__':
    main()

    
