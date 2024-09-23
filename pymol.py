import csv
import argparse
from pymol import cmd, finish_launching

def parse_csv(csv_file):
    patch_rank_0 = []
    all_residues = []
    with open(csv_file, 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            residue_spec = f"{row['residue_id']}/{row['chain']}"  # Correct format for PyMOL selection
            all_residues.append(residue_spec)
            if row['patch_rank'] == '0':
                patch_rank_0.append(residue_spec)
    return patch_rank_0, all_residues

def color_residues(pdb_file, csv_file):
    patch_rank_0, all_residues = parse_csv(csv_file)
    print(all_residues)
    
    cmd.load(pdb_file, "structure")
    cmd.color("blue", "structure")  # Set the initial color for the entire structure
    
    for residue in all_residues:
        try:
            cmd.color("gray", f"chain {residue.split('/')[1]} and resi {residue.split('/')[0]}")
        except:
            print(f"Error coloring residue: {residue}")
    
    for residue in patch_rank_0:
        try:
            cmd.color("red", f"chain {residue.split('/')[1]} and resi {residue.split('/')[0]}")
        except:
            print(f"Error coloring residue: {residue}")
    
    cmd.show("surface", "structure")
    cmd.hide("ribbon", "structure")
    cmd.hide("everything", "structure")
    cmd.show("surface", "structure")

def main():
    parser = argparse.ArgumentParser(description='Color residues in a PDB file based on a CSV file')
    parser.add_argument('pdb_file', type=str, help='Path to the PDB file')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file')

    args = parser.parse_args()
    
    finish_launching()  # Initialize PyMOL
    color_residues(args.pdb_file, args.csv_file)
    cmd.save('output.pse')  # Save the session to an output file if needed

if __name__ == '__main__':
    main()
