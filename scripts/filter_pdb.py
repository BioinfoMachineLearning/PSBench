import os
import shutil
import copy
from util import makedir_if_not_exists
from filter_pdb_util import (
    split_pdb,
    parse_fasta,
    make_chain_id_map,
    get_sequence,
    align_casp_to_native,
    align_and_filter_residues,
    reindex_pdb_files,
    combine_atom_file,
    compute_chain_mapping,
    remove_pdb_file_head,
    PDB_CHAIN_IDS
)

def filter_and_reindex_pdbs(
    fasta: str,
    predictedpdbs: str,
    nativepdb: str,
    outdir: str,
    tmpdir: str,
    clustalw: str
):
    # Parse target name from fasta file
    targetname = os.path.splitext(os.path.basename(fasta))[0]

    makedir_if_not_exists(tmpdir)
    makedir_if_not_exists(outdir)

    pdboutdir = os.path.join(outdir, 'pdb_filtered')
    predoutdir = os.path.join(outdir, 'pred_filtered')
    makedir_if_not_exists(pdboutdir)
    makedir_if_not_exists(predoutdir)

    tmpdir_target = os.path.join(tmpdir, targetname)
    makedir_if_not_exists(tmpdir_target)
    makedir_if_not_exists(os.path.join(tmpdir_target, 'native'))

    # Process native
    true_pdb_file = nativepdb
    true_pdb_chain_files = split_pdb(true_pdb_file, os.path.join(tmpdir_target, 'native'))

    with open(fasta) as f:
        input_fasta_str = f.read()
    input_seqs, input_descs = parse_fasta(input_fasta_str)
    fasta_chains = make_chain_id_map(sequences=input_seqs, descriptions=input_descs)

    native_chains = {cid: get_sequence(f) for cid, f in true_pdb_chain_files.items()}
    fasta_to_native_mapping = align_casp_to_native(fasta_chains, native_chains, clustalw, os.path.join(tmpdir_target, 'native'))

    # print("FASTA to native chain mapping:", fasta_to_native_mapping)

    
    # if os.path.exists(native_filtered_save_dir")
    true_pdb_file_reindex_list = []
    for fasta_cid in fasta_to_native_mapping:
        native_cid = fasta_to_native_mapping[fasta_cid]
        casp_seq = fasta_chains[fasta_cid]
        native_seq_res = native_chains[native_cid]
        pdb_seq_res = copy.deepcopy(native_seq_res)

        native_indices_keep, pdb_indices_keep, _ = align_and_filter_residues(
            clustalw, casp_seq, native_seq_res, pdb_seq_res, os.path.join(tmpdir_target, 'native'), 'native', native_cid
        )

        true_reindexed, _ = reindex_pdb_files(
            true_pdb_chain_files[native_cid], true_pdb_chain_files[native_cid],
            os.path.join(tmpdir_target, 'native'), targetname, native_cid,
            native_indices_keep, pdb_indices_keep
        )
        true_pdb_file_reindex_list.append(true_reindexed)

    true_pdb_final = os.path.join(tmpdir_target, f"{targetname}.pdb")
    combine_atom_file(true_pdb_file_reindex_list, true_pdb_final)
    # native_filtered_save_dir = 
    if not os.path.exists(os.path.join(pdboutdir, f"{targetname}.pdb")):
        shutil.copy(true_pdb_final, os.path.join(pdboutdir, f"{targetname}.pdb"))
    else:
        print(f"File {os.path.join(pdboutdir, f'{targetname}.pdb')} has already been generated")

    for pred_file in os.listdir(predictedpdbs):
        ori_pred_pdb_path = os.path.join(predictedpdbs, pred_file)
        if os.path.isdir(ori_pred_pdb_path):
            continue

        print(f"Processing predicted model: {ori_pred_pdb_path}")
        pred_name = os.path.splitext(pred_file)[0]
        pred_pdb_file = os.path.join(tmpdir_target, f"{pred_name}.pdb")
        remove_pdb_file_head(ori_pred_pdb_path, pred_pdb_file)

        predoutdir_target = os.path.join(predoutdir, targetname)
        makedir_if_not_exists(predoutdir_target)
        if not os.path.exists(os.path.join(predoutdir_target, f"{pred_name}.pdb")):
        # shutil.copy(pred_pdb_final, os.path.join(predoutdir_target, f"{pred_name}.pdb"))

            model_workdir = os.path.join(tmpdir_target, pred_name)
            makedir_if_not_exists(model_workdir)
            pred_pdb_chain_files = split_pdb(pred_pdb_file, model_workdir)

            predicted_chains = {cid: get_sequence(f) for cid, f in pred_pdb_chain_files.items()}
            predicted_to_fasta_mapping = compute_chain_mapping(predicted_chains, fasta_chains, clustalw, model_workdir)

            # print("Predicted to FASTA chain mapping:", predicted_to_fasta_mapping)

            pred_pdb_file_reindex_list = []
            new_chain_ids = iter(PDB_CHAIN_IDS)

            for pred_cid, fasta_cid in predicted_to_fasta_mapping.items():
                native_cid = fasta_to_native_mapping[fasta_cid]
                casp_seq = fasta_chains[fasta_cid]
                native_seq_res = native_chains[native_cid]
                pdb_seq_res = get_sequence(pred_pdb_chain_files[pred_cid])

                native_indices_keep, pdb_indices_keep, _ = align_and_filter_residues(
                    clustalw, casp_seq, native_seq_res, pdb_seq_res, model_workdir, pred_name, native_cid
                )

                _, pred_reindexed = reindex_pdb_files(
                    true_pdb_chain_files[native_cid], pred_pdb_chain_files[pred_cid],
                    model_workdir, targetname, next(new_chain_ids),
                    native_indices_keep, pdb_indices_keep
                )
                pred_pdb_file_reindex_list.append(pred_reindexed)

            pred_pdb_final = os.path.join(model_workdir, f"{pred_name}.pdb")
            combine_atom_file(pred_pdb_file_reindex_list, pred_pdb_final)

            # predoutdir_target = os.path.join(predoutdir, targetname)
            # makedir_if_not_exists(predoutdir_target)
            shutil.copy(pred_pdb_final, os.path.join(predoutdir_target, f"{pred_name}.pdb"))
            print(f"Saved filtered model to {os.path.join(predoutdir_target, f'{pred_name}.pdb')}")
        else:
            print(f"File {os.path.join(predoutdir_target, f'{pred_name}.pdb')} has already been generated")
