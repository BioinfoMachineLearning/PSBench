import os, sys, argparse
from multiprocessing import Pool
from util import is_file, is_dir, makedir_if_not_exists
from filter_pdb import filter_and_reindex_pdbs
from util import extract_tmscore_file
from tqdm import tqdm
import json
import csv


def run_all_commands(params):
    usalign_program, ema_exec, pdb1, pdb2, usalign_outfile, openstructure_outfile = params

    # USAlign command
    usalign_cmd = f"{usalign_program} -ter 1 -tmscore 6 {pdb1} {pdb2} > {usalign_outfile}"
    print("Running USAlign:", usalign_cmd)
    os.system(usalign_cmd)

    # OpenStructure command
    openstruct_cmd = [ema_exec, "compare-structures",
                      "-m", pdb1, "-r", pdb2,
                      "-mf", "pdb", "-rf", "pdb", "-rna", "--out", openstructure_outfile,
                      "--ics", "--ips", "--qs-score", "--lddt", "--local-lddt",
                      "--rigid-scores", "--patch-scores", "--tm-score",
                      "--dockq", "--ics-trimmed", "--ips-trimmed"]
    openstruct_cmd = ' '.join(openstruct_cmd)
    print("Running OpenStructure:", openstruct_cmd)
    os.system(openstruct_cmd)

def run_usalign_only(params):
    usalign_program, pdb1, pdb2, outfile = params
    cmd = f"{usalign_program} -ter 1 -tmscore 6 {pdb1} {pdb2} > {outfile}"
    print("Running USAlign (filtered):", cmd)
    os.system(cmd)



def generate_csv(model_dir,result_dir,output_csv_path):
    content = [["model_name","ics","ics_precision","ics_recall","ips","qs_global","qs_best","lddt","tmscore_mmalign","rmsd","dockq_wave","tmscore_usalign","tmscore_usalign_aligned"]]
    for model in tqdm(os.listdir(model_dir)):
        model_name = model
        openstructure_path = os.path.join(result_dir,f'{model.split(".")[0]}.pdb_openstructure_out')
        tmscore_usalign_path = os.path.join(result_dir,f'{model.split(".")[0]}.pdb_usalign_out')
        tmscore_usalign_aligned_path = os.path.join(result_dir,f'{model.split(".")[0]}.pdb_filt_usalign_out')
        
        # if not os.path.exists(openstructure_path) or not os.path.exists(tmscore_usalign_path) or not os.path.exists(tmscore_usalign_aligned_path):
        #     print("Incomplete data")
        #     continue
        if not is_file(openstructure_path):
            print(f"OpenStructure output file {openstructure_path} does not exist")
            continue
        if not is_file(tmscore_usalign_path):
            print(f"USalign output file {tmscore_usalign_path} does not exist")
            continue
        if not is_file(tmscore_usalign_aligned_path):
            print(f"USalign output file {tmscore_usalign_aligned_path} does not exist")

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


    with open(output_csv_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(content)

    print("saved: ",output_csv_path)
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--fasta', type=is_file, required=True, help="Directory of Fasta file")
    parser.add_argument('--indir', type=is_dir, required=True, help="Directory of original model PDBs")
    parser.add_argument('--nativedir', type=is_file, required=True, help="Native PDB for original models")
    parser.add_argument('--outdir', type=is_dir, required=True, help="Main output directory")
    parser.add_argument('--usalign_program', type=is_file, required=True, help="Path to USAlign executable")
    parser.add_argument('--clustalw_program', type=is_file, required=True, help="Path to USAlign executable")
    parser.add_argument('--nproc', type=int, default=10, help="Number of parallel processes")

    args = parser.parse_args()
    targetname = os.path.basename(os.path.normpath(args.indir))
    native_pdb = args.nativedir
    outdir = os.path.join(args.outdir, targetname)
    output_csv_path = os.path.join(outdir,f"{targetname}.csv")
    makedir_if_not_exists(outdir)

    ema_exec = (
        'docker run --rm -v $(pwd):/home '
        f'-v {args.indir}:{args.indir} '
        f'-v {args.nativedir}:{args.nativedir} '
        f'-v {args.outdir}:{args.outdir} '
        '-u $(id -u $USER):$(id -g $USER) '
        'registry.scicore.unibas.ch/schwede/openstructure:latest'
        )

    # print(ema_exec)

    


    ##run filtration
    tempdir = os.path.join(outdir,"temp")
    makedir_if_not_exists(tempdir)
    filtered_path = os.path.join(outdir,'filtered_pdbs',targetname)
    makedir_if_not_exists(filtered_path)
    filter_and_reindex_pdbs(args.fasta,args.indir,args.nativedir,filtered_path,tempdir,args.clustalw_program)
    indir_filt = os.path.join(filtered_path,"pred_filtered",targetname)
    native_pdb_filt = os.path.join(filtered_path,"pdb_filtered",f"{targetname}.pdb")

    # run OpenStructure and USalign
    outdir = os.path.join(outdir, "results")
    # === Run USAlign + OpenStructure on original ===
    process_list = []
    for model in os.listdir(args.indir):
        infile = os.path.join(args.indir, model)
        usalign_outfile = os.path.join(outdir, f"{model}_usalign_out")
        openstruct_outfile = os.path.join(outdir, f"{model}_openstructure_out")
        if os.path.exists(openstruct_outfile) and os.path.exists(usalign_outfile):
            continue

        process_list.append([args.usalign_program, ema_exec, infile, native_pdb, usalign_outfile, openstruct_outfile])

    #Run only USAlign on filtered
    filt_process_list = []
    if indir_filt and native_pdb_filt:
        targetname_filt = os.path.basename(os.path.normpath(indir_filt))
        outdir_filt = outdir
        makedir_if_not_exists(outdir_filt)

        for model in os.listdir(indir_filt):
            infile = os.path.join(indir_filt, model)
            usalign_outfile = os.path.join(outdir_filt, f"{model}_filt_usalign_out")
            if os.path.exists(usalign_outfile):
                continue
            filt_process_list.append([args.usalign_program, infile, native_pdb_filt, usalign_outfile])

    # Run USAlign + OpenStructure jobs
    pool = Pool(processes=args.nproc)
    pool.map(run_all_commands, process_list)
    pool.close()
    pool.join()

    # Run filtered USAlign jobs
    if filt_process_list:
        pool = Pool(processes=args.nproc)
        pool.map(run_usalign_only, filt_process_list)
        pool.close()
        pool.join()
    
    
    generate_csv(args.indir,outdir,output_csv_path)

