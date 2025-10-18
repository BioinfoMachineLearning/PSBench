
import gemmi
import os
import string
import sys
# import os
from collections import defaultdict


# Standard protein residues (20 amino acids + common modified ones)
PROTEIN_RESIDUES = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    # Modified/non-standard amino acids commonly found
    'MSE', 'SEP', 'TPO', 'PTR', 'HYP', 'MLY', 'CSO', 'CSD', 'CSE',
    'SEC', 'PYL', 'UNK'
}

# DNA residues
DNA_RESIDUES = {'DA', 'DT', 'DG', 'DC', 'DI', 'DU'}

# RNA residues
RNA_RESIDUES = {'A', 'U', 'G', 'C', 'I', 'T'}

# Common ligands and water
EXCLUDE_RESIDUES = {'HOH', 'WAT', 'SO4', 'PO4', 'GOL', 'EDO', 'ACT', 'NAG', 'MAN', 'FUC'}


def filter_protein_only(input_pdb, output_pdb, verbose=True):
    """
    Filter PDB to keep only protein residues.
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
        verbose: Print statistics
    """
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    kept_lines = []
    excluded_residues = set()
    stats = {
        'protein': 0,
        'dna': 0,
        'rna': 0,
        'ligand': 0,
        'water': 0,
        'other': 0
    }
    
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            res_name = line[17:20].strip()
            
            if res_name in PROTEIN_RESIDUES:
                kept_lines.append(line)
                stats['protein'] += 1
            elif res_name in DNA_RESIDUES:
                excluded_residues.add(res_name)
                stats['dna'] += 1
            elif res_name in RNA_RESIDUES:
                excluded_residues.add(res_name)
                stats['rna'] += 1
            elif res_name in {'HOH', 'WAT'}:
                excluded_residues.add(res_name)
                stats['water'] += 1
            elif res_name in EXCLUDE_RESIDUES:
                excluded_residues.add(res_name)
                stats['ligand'] += 1
            else:
                excluded_residues.add(res_name)
                stats['other'] += 1
        
        elif line.startswith(('HEADER', 'TITLE', 'COMPND', 'SOURCE', 'KEYWDS', 
                             'EXPDTA', 'AUTHOR', 'REVDAT', 'JRNL', 'REMARK',
                             'DBREF', 'SEQRES', 'MODRES', 'CRYST1', 'ORIGX',
                             'SCALE', 'MTRIX', 'MODEL', 'ENDMDL', 'TER', 'END')):
            kept_lines.append(line)
    
    with open(output_pdb, 'w') as f:
        f.writelines(kept_lines)
    
    if verbose:
        print(f"\nFiltering: {os.path.basename(input_pdb)}")
        print(f"  Protein atoms kept:  {stats['protein']:>6}")
        if stats['dna'] > 0:
            print(f"  DNA atoms removed:   {stats['dna']:>6}")
        if stats['rna'] > 0:
            print(f"  RNA atoms removed:   {stats['rna']:>6}")
        if stats['water'] > 0:
            print(f"  Water atoms removed: {stats['water']:>6}")
        if stats['ligand'] > 0:
            print(f"  Ligand atoms removed:{stats['ligand']:>6}")
        if stats['other'] > 0:
            print(f"  Other atoms removed: {stats['other']:>6}")
        
        if excluded_residues:
            print(f"  Excluded residues: {', '.join(sorted(excluded_residues))}")
        
        print(f"  Output: {output_pdb}")
    
    return stats


# Usage
# to remove entries that are not proteins: : 
# filter_protein_only("path/to/input_pdb_file","path/to/output_pdb_file")