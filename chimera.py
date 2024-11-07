from chimera import runCommand
import csv
import sys

def parse_csv(csv_file):
    patch_rank_0 = []
    low_plddt = []
    all_residues = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            residue_spec = "{}.{}".format(row['residue_id'], row['chain'])
            all_residues.append(residue_spec)
            print(row['plddt'])
            if row['patch_rank'] == '0' and row['plddt'] and float(row['plddt']) < 70:
                low_plddt.append(residue_spec)
            elif row['patch_rank'] == '0':
                patch_rank_0.append(residue_spec)

    return patch_rank_0, low_plddt, all_residues

def color_residues(pdb_file, csv_file):
    patch_rank_0, low_plddt, all_residues = parse_csv(csv_file)
    print(all_residues)
    runCommand('open {}'.format(pdb_file))
    runCommand('color blue')
    
    for residue in all_residues:
        runCommand('color gray :{}'.format(residue))
    
    for residue in patch_rank_0:
        runCommand('color red :{}'.format(residue))

    for residue in low_plddt:
        runCommand('color yellow :{}'.format(residue))
    
    runCommand('surface')
    runCommand('~ribbon')
    runCommand('~display')
    runCommand('wait')

def main():
    color_residues(sys.argv[1], sys.argv[2])

if __name__ == '__main__':
    main()
