
import gemmi
import os
import string
import sys
# import os
from collections import defaultdict


def cif2pdb_convert(input_cif,output_pdb):
    """
    Convert cif file to pdb file while keeping the residue indices intact
    
    Args:
        input_cif: Input CIF file path
        output_pdb: Output PDB file path
    """
    # Load structure
    st = gemmi.read_structure(input_cif)

    # --- Use label_seq_id for residue numbering ---
    for model in st:
        for chain in model:
            for res in chain:
                if res.label_seq is not None:
                    res.seqid.num = res.label_seq

    # Write temporary full PDB
    temp_output = output_pdb + ".tmp"
    st.write_pdb(temp_output)

    # Filter only ATOM, TER, and END lines
    with open(temp_output, "r") as fin, open(output_pdb, "w") as fout:
        for line in fin:
            if line.startswith(("ATOM", "TER", "END")):
                fout.write(line)

    # Remove temporary file
    import os
    os.remove(temp_output)

    print(f"✅ Converted {input_cif} → {output_pdb}")
    print("   - label_seq_id numbering preserved")
    print("   - only ATOM/TER/END lines kept")










# Usage
# to convert a cif file to pdb file: 
# cif2pdb_convert("path/to/input_cif_file","path/to/output_pdb_file")






