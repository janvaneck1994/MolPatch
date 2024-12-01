from chimera import runCommand
import csv
import sys

def parse_csv(csv_file):
    """
    Parses the CSV file to extract residue information.
    """
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

def parse_pdb_b_factors(pdb_file):
    """
    Parses the PDB file to extract residues with low B-factors and matching allowed residue names.
    """
    low_b_factor = []
    allowed_residues = ["ALA", "CYS", "PHE", "ILE", "LEU", "MET", "VAL", "TRP", "TYR"]

    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM"):
                residue_name = line[17:20].strip()  # Extract the residue name (3-letter code)
                chain = line[21].strip()
                residue_id = line[22:26].strip()
                b_factor = float(line[60:66].strip())
                residue_spec = "{}.{}".format(residue_id, chain)

                if (
                    residue_name in allowed_residues  # Check if the residue is in the allowed list
                    and b_factor < 70  # Check if the B-factor is less than 70
                    and residue_spec not in low_b_factor  # Avoid duplicates
                ):
                    low_b_factor.append(residue_spec)

    return low_b_factor


def color_residues(pdb_file, csv_file):
    """
    Colors residues based on information from the CSV and PDB files.
    """
    patch_rank_0, all_residues = parse_csv(csv_file)
    low_b_factor = parse_pdb_b_factors(pdb_file)
    
    runCommand('open {}'.format(pdb_file))
    runCommand('color blue')

    # Color all residues gray initially
    for residue in all_residues:
        runCommand('color gray :{}'.format(residue))

    # Color residues from CSV
    for residue in patch_rank_0:
        runCommand('color red :{}'.format(residue))

    # Color residues with low B-factors from PDB
    for residue in low_b_factor:
        runCommand('color yellow :{}'.format(residue))
    
    runCommand('surface')
    runCommand('~ribbon')
    runCommand('~display')
    runCommand('wait')

def main():
    """
    Main function to run the script.
    """
    pdb_file = sys.argv[1]
    csv_file = sys.argv[2]
    color_residues(pdb_file, csv_file)

if __name__ == '__main__':
    main()
