import chimera
from chimera import runCommand
import csv
import argparse

def parse_csv(csv_file):
    patch_rank_0 = []
    all_residues = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            residue_spec = "{}.{}".format(row['residue_id'], row['chain'])
            all_residues.append(residue_spec)
            if row['patch_rank'] == '0':
                patch_rank_0.append(residue_spec)
    return patch_rank_0, all_residues

def color_residues(pdb_file, csv_file):
    patch_rank_0, all_residues = parse_csv(csv_file)
    print(all_residues)
    runCommand('open {}'.format(pdb_file))
    runCommand('color blue')
    
    for residue in all_residues:
        runCommand('color gray :{}'.format(residue))
    
    for residue in patch_rank_0:
        runCommand('color red :{}'.format(residue))
    
    runCommand('surface')
    runCommand('~ribbon')
    runCommand('~display')
    runCommand('wait')

def main():
    parser = argparse.ArgumentParser(description='Color residues in a PDB file based on a CSV file')
    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file')

    args = parser.parse_args()
    
    color_residues(args.pdb_file, args.csv_file)

if __name__ == '__main__':
    main()
