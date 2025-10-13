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


import os, sys
import numpy as np
import argparse
from pdb_process import *
from protein_util import read_qa_txt_as_df, parse_fasta, complete_result
from util import is_file, is_dir, makedir_if_not_exists, check_contents, read_option_file, check_dirs
from itertools import chain
import shutil
import string
import copy
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

def make_chain_id_map(sequences, descriptions):
    """
    Generate a mapping of chain IDs to sequences from parsed FASTA data.
    
    Args:
        sequences (list): A list of sequences from the FASTA file.
        descriptions (list): A list of descriptions corresponding to each sequence.
        
    Returns:
        dict: A dictionary mapping chain IDs to sequences.
    """
    chain_id_map = {}
    for desc, seq, chain_id in zip(descriptions, sequences, PDB_CHAIN_IDS):
        # Extract the chain ID from the description, assuming the description
        # format is something like ">chain_id description" or ">chain_id"
        chain_id_map[chain_id] = seq
    return chain_id_map

#starting_resid = -1
def renumber_residues(fhandle, starting_resid):
    """Resets the residue number column to start from a specific number.
    """
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    resid = starting_resid - 1  # account for first residue
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                resid += 1
                if resid > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(resid).rjust(4) + line[26:]

        else:
            yield line

#starting_resid = -1
def renumber_residues_with_order(fhandle, reorder_indices):
    """Resets the residue number column to start from a specific number.
    """
    _pad_line = pad_line
    prev_resid = None  # tracks chain and resid
    counter = -1
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        line = _pad_line(line)
        if line.startswith(records):
            line_resuid = line[17:27]
            if line_resuid != prev_resid:
                prev_resid = line_resuid
                counter += 1
                if reorder_indices[counter] > 9999:
                    emsg = 'Cannot set residue number above 9999.\n'
                    sys.stderr.write(emsg)
                    sys.exit(1)

            yield line[:22] + str(reorder_indices[counter]).rjust(4) + line[26:]

        else:
            yield line


#starting_value = -1
def renumber_atom_serials(fhandle, starting_value):
    """Resets the atom serial number column to start from a specific number.
    """

    # CONECT 1179  746 1184 1195 1203
    fmt_CONECT = "CONECT{:>5s}{:>5s}{:>5s}{:>5s}{:>5s}" + " " * 49 + "\n"
    char_ranges = (slice(6, 11), slice(11, 16),
                   slice(16, 21), slice(21, 26), slice(26, 31))

    serial_equiv = {'': ''}  # store for conect statements

    serial = starting_value
    records = ('ATOM', 'HETATM')
    for line in fhandle:
        if line.startswith(records):
            serial_equiv[line[6:11].strip()] = serial
            yield line[:6] + str(serial).rjust(5) + line[11:]
            serial += 1
            if serial > 99999:
                emsg = 'Cannot set atom serial number above 99999.\n'
                sys.stderr.write(emsg)
                sys.exit(1)

        elif line.startswith('ANISOU'):
            # Keep atom id as previous atom
            yield line[:6] + str(serial - 1).rjust(5) + line[11:]

        elif line.startswith('CONECT'):
            # 6:11, 11:16, 16:21, 21:26, 26:31
            serials = [line[cr].strip() for cr in char_ranges]

            # If not found, return default
            new_serials = [str(serial_equiv.get(s, s)) for s in serials]
            conect_line = fmt_CONECT.format(*new_serials)

            yield conect_line
            continue

        elif line.startswith('MODEL'):
            serial = starting_value
            yield line

        elif line.startswith('TER'):
            yield line[:6] + str(serial).rjust(5) + line[11:]
            serial += 1

        else:
            yield line


def keep_residues(fhandle, residue_range):
    """Deletes residues within a certain numbering range.
    """
    prev_res = None
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    for line in fhandle:
        if line.startswith(records):

            res_id = line[21:26]  # include chain ID
            if res_id != prev_res:
                prev_res = res_id

            if int(line[22:26]) not in residue_range:
                continue

        yield line


def reindex_pdb_file(in_file, out_file, keep_indices, reorder_indices = []):
    # print(keep_indices)
    fhandle = open(in_file, 'r')
    # fhandle = renumber_residues(fhandle, 1)
    # fhandle = renumber_atom_serials(fhandle, 1)
    fhandle = keep_residues(fhandle, keep_indices)
    if len(reorder_indices) > 0:
        fhandle = renumber_residues_with_order(fhandle, reorder_indices)
    else:
        fhandle = renumber_residues(fhandle, 1)
    fhandle = renumber_atom_serials(fhandle, 1)
    write_pdb_file(fhandle, out_file)


# s input str
# c search char
def find_all(s, c):
    idx = s.find(c)
    while idx != -1:
        yield idx
        idx = s.find(c, idx + 1)

def get_sequence(inpdb):
    """Enclosing logic in a function to simplify code"""

    seq_to_res_mapping = []
    res_codes = [
        # 20 canonical amino acids
        ('CYS', 'C'), ('ASP', 'D'), ('SER', 'S'), ('GLN', 'Q'),
        ('LYS', 'K'), ('ILE', 'I'), ('PRO', 'P'), ('THR', 'T'),
        ('PHE', 'F'), ('ASN', 'N'), ('GLY', 'G'), ('HIS', 'H'),
        ('LEU', 'L'), ('ARG', 'R'), ('TRP', 'W'), ('ALA', 'A'),
        ('VAL', 'V'), ('GLU', 'E'), ('TYR', 'Y'), ('MET', 'M'),
        # Non-canonical amino acids
        # ('MSE', 'M'), ('SOC', 'C'),
        # Canonical xNA
        ('  U', 'U'), ('  A', 'A'), ('  G', 'G'), ('  C', 'C'),
        ('  T', 'T'),
    ]

    three_to_one = dict(res_codes)
    # _records = set(['ATOM  ', 'HETATM'])
    _records = set(['ATOM  '])

    sequence = []
    read = set()
    for line in open(inpdb):
        line = line.strip()
        if line[0:6] in _records:
            resn = line[17:20]
            resi = line[22:26]
            icode = line[26]
            r_uid = (resn, resi, icode)
            if r_uid not in read:
                read.add(r_uid)
            else:
                continue
            aa_resn = three_to_one.get(resn, 'X')
            sequence.append(aa_resn)
            seq_to_res_mapping += [int(resi)]

    return {'sequence': ''.join(sequence), 'mapping': seq_to_res_mapping}


def split_pdb(complex_pdb, outdir):
    pdbs = {}
    pre_chain = None
    i = 0
    for line in open(complex_pdb, 'r').readlines():
        if not line.startswith('ATOM'):
            continue
        chain_name = line[21]
        if pre_chain is None:
            pre_chain = chain_name
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
        elif chain_name == pre_chain:
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
        else:
            fw.close()
            i = i + 1
            fw = open(outdir + '/' + chain_name + '.pdb', 'w')
            pdbs[chain_name] = outdir + '/' + chain_name + '.pdb'
            # fw.write(line[:21] + ' ' + line[22:])
            fw.write(line)
            pre_chain = chain_name
    fw.close()
    return pdbs

def remove_head(fhandle):
    """Remove head information.
    """
    records = ('ATOM', 'TER', 'END')
    for line in fhandle:
        if not line.startswith(records):
            continue
        yield line

def remove_pdb_file_head(in_file, out_file):
    fhandle = open(in_file, 'r')
    fhandle = remove_head(fhandle)
    write_pdb_file(fhandle, out_file)


def write_pdb_file(new_pdb, pdb_file):
    if os.path.exists(pdb_file):
        os.remove(pdb_file)
    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                open(pdb_file, 'a').write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        open(pdb_file, 'a').write(''.join(_buffer))
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

def align_sequences(clustalw_program, casp_seq, native_seq, pdb_seq, outfile, workdir):

    with open(workdir + '/' + outfile, 'w') as fw:
        if len(pdb_seq) == 0:
            fw.write(f"%CASP\n{casp_seq}\n%NATIVE\n{native_seq}")
        else:    
            fw.write(f"%CASP\n{casp_seq}\n%NATIVE\n{native_seq}\n%PDB\n{pdb_seq}")

    cmd = f"{clustalw_program}  -MATRIX=BLOSUM -TYPE=PROTEIN -INFILE={workdir}/{outfile} -OUTFILE={workdir}/{outfile}.out >/dev/null 2>&1"
    os.system(cmd)

    casp_align_seq = ""
    native_align_seq = ""
    pdb_align_seq = ""
    for line in open(workdir + '/' + outfile + '.out'):
        if line[0:4] == "CASP":
            casp_align_seq += line.split()[1].rstrip('\n')
        if line[0:6] == "NATIVE":
            native_align_seq += line.split()[1].rstrip('\n')
        if line[0:3] == "PDB":
            pdb_align_seq += line.split()[1].rstrip('\n')

    return casp_align_seq, native_align_seq, pdb_align_seq

def cal_sequence_identity_from_seq(seq1, seq2):
    common_count = 0
    for i in range(len(seq1)):
        char1 = seq1[i:i + 1]
        char2 = seq2[i:i + 1]
        if char1 != '-' and char2 != 'X' and char1 == char2:
            common_count += 1
    seqid = float(common_count) / float(len([1 for i in range(len(seq2)) if seq2[i] != '-' and seq2[i] != 'X']))
    return seqid


def cal_sequence_identity(clustalw_program, seq1, seq2, workdir):
    with open(workdir + '/tmp.fasta', 'w') as fw:
        fw.write(f"%SEQ1\n{seq1}\n%SEQ2\n{seq2}")

    cmd = f"{clustalw_program}  -MATRIX=BLOSUM -TYPE=PROTEIN -INFILE={workdir}/tmp.fasta -OUTFILE={workdir}/tmp.fasta.out >/dev/null 2>&1"
    os.system(cmd)
    # print(cmd)

    seq1 = ""
    seq2 = ""
    for line in open(workdir + '/tmp.fasta.out'):
        if line[0:4] == "SEQ1":
            seq1 += line.split()[1].rstrip('\n')
        if line[0:4] == "SEQ2":
            seq2 += line.split()[1].rstrip('\n')

    return cal_sequence_identity_from_seq(seq1, seq2)

def combine_atom_file(atom_file_list, complex_file):
    """Combine multiple atom files into a single PDB file and add 'TER' at the end of each chain."""
    atom = []
    for file in atom_file_list:
        with open(file, 'r') as f:
            lines = f.readlines()
            atom.extend(lines)
            if not lines[-1].startswith("TER"):  # Ensure 'TER' is added if not present
                atom.append("TER\n")  # Add 'TER' to mark the end of the chain

    # Add 'END' to mark the end of the PDB file
    atom.append("END\n")

    # Write the combined content to the complex file
    with open(complex_file, 'w') as myfile:
        myfile.write(''.join(atom))

    print(f"Combined PDB file saved as {complex_file}")

def compute_chain_mapping(predicted_chains, fasta_chains, clustalw_program, workdir):
    """Compute the best one-to-one mapping of predicted chains to CASP FASTA chains based on sequence identity."""
    chain_mapping = {}
    mapped_fasta_chains = set()  # Keep track of already mapped CASP chains

    for pred_chain_id, pred_seq in predicted_chains.items():
        best_match = None
        best_seqid = 0
        for fasta_chain_id, fasta_seq in fasta_chains.items():
            # Calculate sequence identity using the sequence directly
            seqid = cal_sequence_identity(clustalw_program, pred_seq['sequence'], fasta_seq, workdir)
            if seqid > best_seqid:  # Find the best match based on sequence identity
                best_seqid = seqid
                best_match = fasta_chain_id

        if best_match:
            chain_mapping[pred_chain_id] = best_match

    return chain_mapping

def align_casp_to_native(fasta_chains, native_chains, clustalw_program, workdir):
    """Compute the best one-to-one mapping of CASP FASTA chains to native chains based on sequence identity."""
    casp_to_native_mapping = {}
    mapped_native_chains = set()  # Keep track of already mapped native chains

    for fasta_chain_id, fasta_seq in fasta_chains.items():
        best_match = None
        best_seqid = 0
        for native_chain_id, native_seq in native_chains.items():
            if native_chain_id in mapped_native_chains:
                continue  # Skip if this native chain is already mapped

            # Calculate sequence identity
            seqid = cal_sequence_identity(clustalw_program, fasta_seq, native_seq['sequence'], workdir)
            if seqid > best_seqid:  # Find the best match based on sequence identity
                best_seqid = seqid
                best_match = native_chain_id

        if best_match:
            casp_to_native_mapping[fasta_chain_id] = best_match
            mapped_native_chains.add(best_match)  # Mark this native chain as mapped

    return casp_to_native_mapping

def reorder_predicted_chains(pred_pdb_chain_files, predicted_to_fasta_mapping, fasta_to_native_mapping):
    """
    Reorder predicted PDB chains based on the combined mapping through the CASP FASTA file.

    Args:
        pred_pdb_chain_files (dict): Dictionary of predicted PDB chain files.
        predicted_to_fasta_mapping (dict): Mapping of predicted chains to CASP FASTA chains.
        fasta_to_native_mapping (dict): Mapping of CASP FASTA chains to native chains.

    Returns:
        dict: Reordered dictionary where keys are native chain IDs and values are predicted PDB chain files.
    """
    reordered_chains = {}

    for pred_chain_id, fasta_chain_id in predicted_to_fasta_mapping.items():
        # Find the native chain ID corresponding to the mapped FASTA chain
        native_chain_id = fasta_to_native_mapping.get(fasta_chain_id)
        
        if native_chain_id:
            # Ensure the predicted chain exists in the provided chain files
            if pred_chain_id in pred_pdb_chain_files:
                reordered_chains[native_chain_id] = pred_pdb_chain_files[pred_chain_id]
            else:
                print(f"Warning: Predicted chain {pred_chain_id} not found in predicted PDB files.")
        else:
            print(f"Warning: No mapping found for FASTA chain {fasta_chain_id} in native chains.")

    return reordered_chains

def align_and_filter_residues(clustalw_program, casp_seq, native_seq_res, pdb_seq_res, model_workdir, pred_name, chain_id):
    """
    Align sequences and determine the residues to keep for both the native and predicted PDB files.
    
    Args:
        clustalw_program (str): Path to the CLUSTALW program for sequence alignment.
        casp_seq (str): CASP sequence from the FASTA file.
        native_seq_res (dict): Native sequence and residue mapping.
        pdb_seq_res (dict): Predicted PDB sequence and residue mapping.
        model_workdir (str): Working directory for the model.
        pred_name (str): Name of the prediction.
        chain_id (str): Chain identifier.
        
    Returns:
        tuple: Lists of native indices to keep, PDB indices to keep, and reordered PDB indices.
    """
    # Perform sequence alignment
    casp_align_seq, native_align_seq, pdb_align_seq = align_sequences(
        clustalw_program, casp_seq, native_seq_res['sequence'], pdb_seq_res['sequence'], 
        f"{pred_name}_{chain_id}.fasta", model_workdir
    )
    # print(casp_align_seq)
    # print(native_align_seq)
    # print(pdb_align_seq)
    # Identify positions with gaps in CASP or native alignments
    index_casp_temp = list(find_all(casp_align_seq, '-'))
    index_native_temp = list(find_all(native_align_seq, '-'))
    
    index_not_match = [i for i in range(len(casp_align_seq)) if casp_align_seq[i] != '-' and native_align_seq[i] != '-' and casp_align_seq[i] != native_align_seq[i]]
    # Determine positions to keep (those without gaps in either sequence)
    indexDel = index_casp_temp + index_native_temp + index_not_match
    indices_keep = [i for i in range(len(native_align_seq)) if i not in indexDel]

    # Map indices for native residues to keep
    valid_char_counter = -1
    native_indices_keep = []
    for i in range(len(casp_align_seq)):
        if native_align_seq[i] != '-':
            valid_char_counter += 1
        if i in indices_keep:
            native_indices_keep.append(native_seq_res['mapping'][valid_char_counter])

    # Map indices for predicted PDB residues to keep and reorder
    pdb_indices_keep = []
    pdb_indices_order = []
    pdb_indices_counter = 0
    valid_char_counter = -1
    for i in range(len(casp_align_seq)):
        if native_align_seq[i] != '-':
            pdb_indices_counter += 1
        if pdb_align_seq[i] != '-':
            valid_char_counter += 1
        if i in indices_keep and pdb_align_seq[i] != '-' and native_align_seq[i] != '-':
            pdb_indices_keep.append(pdb_seq_res['mapping'][valid_char_counter])
            pdb_indices_order.append(pdb_indices_counter)

    return native_indices_keep, pdb_indices_keep, pdb_indices_order

def reindex_pdb_file_with_mapping(in_file, out_file, residue_mapping):
    """
    Reindex a PDB file using a mapping from predicted residue indices to native residue indices.
    
    Args:
        in_file (str): Input PDB file.
        out_file (str): Output PDB file after reindexing.
        residue_mapping (dict): Mapping of original residue indices to new residue indices.
    """
    fhandle = open(in_file, 'r')
    records = ('ATOM', 'HETATM', 'TER', 'ANISOU')
    
    with open(out_file, 'w') as out_f:
        for line in fhandle:
            if line.startswith(records):
                orig_residue_index = int(line[22:26].strip())
                if orig_residue_index in residue_mapping:
                    new_residue_index = residue_mapping[orig_residue_index]
                    line = line[:22] + str(new_residue_index).rjust(4) + line[26:]
                    out_f.write(line)
            else:
                out_f.write(line)

def reindex_pdb_files(true_pdb_file, pred_pdb_file, model_workdir, targetname, chain_id, 
                      native_indices_keep, pdb_indices_keep):
    """
    Reindex the PDB files based on the filtered residues using the native residue indices.
    
    Args:
        true_pdb_file (str): Path to the native PDB file.
        pred_pdb_file (str): Path to the predicted PDB file.
        model_workdir (str): Working directory for the model.
        targetname (str): Name of the target.
        chain_id (str): Chain identifier.
        native_indices_keep (list): List of native residue indices to keep.
        pdb_indices_keep (list): List of predicted PDB residue indices to keep.
        
    Returns:
        tuple: Paths to the reindexed native and predicted PDB files.
    """
    # Reindex the native PDB file based on the native indices to keep
    true_pdb_file_reindex_temp = f'{model_workdir}/{targetname}_{chain_id}_filtered.pdb'
    reindex_pdb_file(true_pdb_file, true_pdb_file_reindex_temp, native_indices_keep, native_indices_keep)
    
    # Map the predicted residue indices to the native residue indices
    # This step aligns the predicted indices with the corresponding native indices
    residue_mapping = {pred_idx: native_idx for pred_idx, native_idx in zip(pdb_indices_keep, native_indices_keep)}
    # print(residue_mapping) 
    # Reindex the predicted PDB file to follow the residue numbering from the native PDB
    pred_pdb_file_reindex_temp = f'{model_workdir}/{chain_id}_filtered.pdb'
    reindex_pdb_file_with_mapping(pred_pdb_file, pred_pdb_file_reindex_temp, residue_mapping)

    return true_pdb_file_reindex_temp, pred_pdb_file_reindex_temp

def reorder_chains_in_pdb(in_file, out_file, chain_order):
    """
    Reorder the chains in a PDB file to match a specified chain order.
    
    Args:
        in_file (str): Path to the input PDB file.
        out_file (str): Path to the output PDB file with reordered chains.
        chain_order (list): List of chain identifiers specifying the desired order.
    """
    # Dictionary to store lines for each chain
    chain_lines = {chain_id: [] for chain_id in chain_order}

    # Read the input PDB file and categorize lines by chain
    with open(in_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM', 'ANISOU', 'TER')):
                chain_id = line[21]  # Extract chain identifier
                if chain_id in chain_lines:
                    chain_lines[chain_id].append(line)
            else:
                # Keep non-atom lines like END
                chain_lines.setdefault('other', []).append(line)

    # Write the reordered chains to the output file
    with open(out_file, 'w') as f:
        for chain_id in chain_order:
            if chain_id in chain_lines:
                f.writelines(chain_lines[chain_id])
        # Write other lines like END at the end
        if 'other' in chain_lines:
            f.writelines(chain_lines['other'])


def reorder_predicted_pdb_chains(native_pdb_file, pred_pdb_file, model_workdir, targetname, pred_name):
    """
    Reorder the chains in the predicted PDB file to match the chain order in the native PDB file.
    
    Args:
        native_pdb_file (str): Path to the native PDB file.
        pred_pdb_file (str): Path to the predicted PDB file.
        model_workdir (str): Working directory for the model.
        targetname (str): Name of the target.
        pred_name (str): Name of the predicted structure.
        
    Returns:
        str: Path to the reordered predicted PDB file.
    """
    # Determine the chain order from the native PDB file
    chain_order = []
    with open(native_pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]  # Extract chain identifier
                if chain_id not in chain_order:
                    chain_order.append(chain_id)

    # Reorder the predicted PDB chains to match the native PDB chain order
    reordered_pred_pdb_file = f'{model_workdir}/{pred_name}_reordered.pdb'
    reorder_chains_in_pdb(pred_pdb_file, reordered_pred_pdb_file, chain_order)

    return reordered_pred_pdb_file