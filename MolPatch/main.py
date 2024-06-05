#!/usr/bin/env python
import argparse
from ResidueBased.residue_based_molpatch import calculate_patches
from preprocessing import preprocess
import os

def main():
    parser = argparse.ArgumentParser(description="Process input file for residue-based molpatch")

    parser.add_argument("-rh", "--remove-hoh", action="store_true", help="Remove water (HOH) molecules")
    parser.add_argument("-rhetm", "--remove_hetm", action="store_true", help="Remove HETM molecules")
    parser.add_argument("-sc", "--split_chains", action="store_true", help="Split chains into separate files")
    parser.add_argument("-alpha", "--alphafold", default=False, action="store_true", help="Input is Alphafold prediction")
    parser.add_argument("-plt", "--plddt_threshold", default=65, type=float, help="The local confidence threshold")
    parser.add_argument("-mpl", "--minimum_plddt", default=5, type=int, help="Minimum residues with a lower plDDT than the set threshold")


    args = parser.parse_args()

    remove_hoh = args.remove_hoh
    remove_hetm = args.remove_hetm
    split_chains = args.split_chains
    plddt_threshold = args.plddt_threshold
    minimum_plddt = args.minimum_plddt
    alphafold = args.alphafold

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
        preprocessed_pdbs = []

        try:
            preprocessed_pdbs = preprocess(infile=input_file, 
                                        output_dir=output_dir, 
                                        split_chains=split_chains, 
                                        remove_hetatm=remove_hetm, 
                                        remove_hoh=remove_hoh,
                                        plddt_threshold=plddt_threshold,
                                        minimum_plddt=minimum_plddt,
                                        alphafold=alphafold
                                        )
        except:
            print(f"Unable to preprocess {filename}")

        
        for file in preprocessed_pdbs:
            try:
                calculate_patches(infile=file, residues=residues, plot=False, output_dir=output_dir)
            except:
                print(f"Unable to calculate {file}")

if __name__ == "__main__":
    main()
