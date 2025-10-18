
import gemmi
import os
import string
import sys
# import os
from collections import defaultdict

def normalize_insertion_codes(input_pdb, output_pdb):
    """
    Remove insertion codes from PDB file and renumber residues sequentially.
    
    Args:
        input_pdb: Input PDB file path
        output_pdb: Output PDB file path
    """
    
    with open(input_pdb, 'r') as f:
        lines = f.readlines()
    
    # Track new residue numbers for each chain
    chain_residue_map = defaultdict(dict)
    chain_last_resnum = defaultdict(int)
    
    # First pass: build mapping of old to new residue numbers
    print("Building residue number mapping...")
    for line in lines:
        if line.startswith(('ATOM', 'HETATM')):
            chain_id = line[21]
            res_num = int(line[22:26].strip())
            ins_code = line[26].strip()
            
            old_key = (res_num, ins_code)
            
            if old_key not in chain_residue_map[chain_id]:
                if ins_code:
                    # Has insertion code - assign next sequential number after base residue
                    chain_last_resnum[chain_id] += 1
                    new_resnum = chain_last_resnum[chain_id]
                    chain_residue_map[chain_id][old_key] = new_resnum
                    print(f"  Chain {chain_id}: {res_num}{ins_code} -> {new_resnum}")
                else:
                    # No insertion code - keep original number
                    chain_residue_map[chain_id][old_key] = res_num
                    chain_last_resnum[chain_id] = res_num
    
    # Second pass: write new PDB with updated residue numbers
    print(f"\nWriting output to {output_pdb}...")
    with open(output_pdb, 'w') as f:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                chain_id = line[21]
                res_num = int(line[22:26].strip())
                ins_code = line[26].strip()
                
                old_key = (res_num, ins_code)
                new_res_num = chain_residue_map[chain_id][old_key]
                
                # Reconstruct line with new residue number and no insertion code
                new_line = (
                    line[:22] +
                    f"{new_res_num:4d}" +
                    " " +
                    line[27:]
                )
                f.write(new_line)
            else:
                # Keep all other lines unchanged (headers, etc.)
                f.write(line)
    
    # Print summary
    total_residues = sum(len(mapping) for mapping in chain_residue_map.values())
    chains_with_insertions = sum(1 for chain, mapping in chain_residue_map.items() 
                                  if any(ins_code for (_, ins_code) in mapping.keys()))
    
    print(f"\nSummary:")
    print(f"  Total chains processed: {len(chain_residue_map)}")
    print(f"  Chains with insertion codes: {chains_with_insertions}")
    print(f"  Total residues: {total_residues}")
    print(f"  Output saved to: {output_pdb}")

# Usage
# to normalize the residue indices to exclude insertion codes : 
# normalize_insertion_codes("path/to/input_pdb_file","path/to/output_pdb_file")