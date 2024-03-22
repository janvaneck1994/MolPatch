#!/usr/bin/env python
import argparse
from ResidueBased.residue_based_molpatch import calculate_patches
from preprocessing import preprocess
import os

def main():
    parser = argparse.ArgumentParser(description="Process input file for residue-based molpatch")

    parser.add_argument("-rh", "--remove-hoh", action="store_true", help="Remove water (HOH) molecules")
    parser.add_argument("-rhetm", "--remove_hetm", action="store_true", help="Remove HETM molecules")
    parser.add_argument("-sc", "--split-chains", action="store_true", help="Split chains into separate files")

    args = parser.parse_args()

    remove_hoh = args.remove_hoh
    remove_hetm = args.remove_hetm
    split_chains = args.split_chains

    residues = ["A","C","F","I","L","M","V","W","Y"]

    input_dir = "/input/"

    for input in os.listdir(input_dir):
        if not input.endswith('.pdb'):
            print(f"File is not a PDB file: {input}")
            continue

        input_file = os.path.join(input_dir, input)
        filename = os.path.basename(input_file)
        pdb_id = filename.split('.')[0]

        output_dir = "/output/" + pdb_id
        os.makedirs(output_dir, exist_ok=True)
        
        preprocessed_pdbs = preprocess(infile=input_file, 
                                    output_dir=output_dir, 
                                    split_chains=split_chains, 
                                    remove_hetatm=remove_hetm, 
                                    remove_hoh=remove_hoh
                                    )

        for file in preprocessed_pdbs:
            calculate_patches(infile=file, residues=residues, plot=False, output_dir=output_dir)

if __name__ == "__main__":
    main()
