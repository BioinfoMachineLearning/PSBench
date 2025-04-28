import os
import json
import csv
from tqdm import tqdm
import argparse
from util import is_dir,is_file
# from Bio import PDB

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

# for target in target_list:
if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-pp','--model_dir', type=is_dir, required=True, help="Directory containing models for one target")
    parser.add_argument('-os','--openstructure_result_dir', type=is_dir, required=True, help="Directory where openstructure results for the target is located")
    parser.add_argument('-tm_u','--tmscore_usalign_results_dir', type=is_dir, required=True, help="Directory where tmscore_usalign results for the target is located")
    parser.add_argument('-tm_ua','--tmscore_usalign_aligned_results_dir', type=is_dir, required=True, help="Directory where tmscore_usalign_aligned results for the target is located")
    parser.add_argument('-oc','--output_csv', required=True, help="Path where output csv is to be saved")
    args = parser.parse_args()

    model_dir = args.model_dir
    openstructure_result_dir = args.openstructure_result_dir
    tmscore_usalign_results_dir = args.tmscore_usalign_results_dir
    tmscore_usalign_aligned_results_dir = args.tmscore_usalign_aligned_results_dir
    output_csv = args.output_csv
    
    # add content header
    content = [["model_name","ics","ics_precision","ics_recall","ips","qs_global","qs_best","lddt","tmscore_mmalign","rmsd","dockq_wave","tmscore_usalign","tmscore_usalign_aligned"]]
    for model in tqdm(os.listdir(model_dir)):
        model_name = model.split("_filtered.pdb")[0]
        openstructure_path = os.path.join(openstructure_result_dir,f'{model.split(".")[0]}.pdb_filtered_out')
        tmscore_usalign_path = os.path.join(tmscore_usalign_results_dir,f'{model.split(".")[0]}.pdb_filtered_out')
        tmscore_usalign_aligned_path = os.path.join(tmscore_usalign_aligned_results_dir,f'{model.split(".")[0]}_filtered.pdb_filtered_out')
        
        if not os.path.exists(openstructure_path) or not os.path.exists(tmscore_usalign_path) or not os.path.exists(tmscore_usalign_aligned_path):
            print("Incomplete data")
            continue

        with open(openstructure_path,"r") as f:
            openstructure_data = json.load(f)

        if openstructure_data["status"]=="FAILURE":
            print("FAILURE in openstructure for file: ",openstructure_path)
            continue

        ics = openstructure_data.get("ics_trimmed", "NA")
        ics_precision = openstructure_data.get("ics_precision_trimmed", "NA")
        ics_recall = openstructure_data.get("ics_recall_trimmed", "NA")
        ips = openstructure_data.get("ips_trimmed", "NA")
        qs_global = openstructure_data.get("qs_global", "NA")
        qs_best = openstructure_data.get("qs_best", "NA")
        lddt = openstructure_data.get("lddt", "NA")
        rmsd = openstructure_data.get("rmsd", "NA")
        dockq_wave = openstructure_data.get("dockq_wave", "NA")

        tmscore_mmalign = openstructure_data.get("tm_score", "NA")
        tmscore_usalign = extract_tmscore_file(tmscore_usalign_path)
        tmscore_usalign_aligned = extract_tmscore_file(tmscore_usalign_aligned_path)

        content.append([model_name,ics,ics_precision,ics_recall,ips,qs_global,qs_best,lddt,tmscore_mmalign,rmsd,dockq_wave,tmscore_usalign,tmscore_usalign_aligned])


    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(content)

    print("saved: ",output_csv)

