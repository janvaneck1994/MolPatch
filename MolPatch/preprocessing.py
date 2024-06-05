#!/usr/bin/python3

import shutil
import tempfile
import os
from Bio.PDB import Select, PDBParser, MMCIFParser, PDBIO
import os

class ChainSelect(Select):
    def __init__(self, chain):
        self.chain = chain

    def accept_chain(self, chain):
        if chain.get_id() == self.chain:
            return 1
        else:
            return 0
        
class ResidueFilterSelect(Select):
    def __init__(self, residues_to_remove):
        self.residues_to_remove = residues_to_remove

    def accept_residue(self, residue):
        if residue in self.residues_to_remove:
            return 0
        else:
            return 1
        

def preprocess_remove_hetatm(input_file):
    """ Remove heteroatoms from PDB file and move the new file to a temporary folder """

    temp_dir = tempfile.mkdtemp()
    filename = os.path.basename(input_file)
    output_file = os.path.join(temp_dir, filename)

    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                if not line.startswith('HETATM'):
                    outfile.write(line)

    return output_file

def preprocess_remove_hoh(input_file):
    """ Remove water molecules from PDB file and move the new file to a temporary folder """

    temp_dir = tempfile.mkdtemp()
    filename = os.path.basename(input_file)
    output_file = os.path.join(temp_dir, filename)

    with open(input_file, 'r') as infile:
        with open(output_file, 'w') as outfile:
            for line in infile:
                if input_file.endswith('.pdb') and line[17:20] != 'HOH':
                    outfile.write(line)
                elif input_file.endswith('.cif') and line[22:25] != 'HOH':
                    outfile.write(line)

    return output_file


def preprocess_split_chains(input_file):
    """ Split PDB file into separate PDB files, each containing one chain of the PDB structure """

    temp_dir = tempfile.mkdtemp()
    pdb_id = os.path.splitext(os.path.basename(input_file))[0]

    if input_file.endswith('.pdb'):
        parser = PDBParser(QUIET=1)
        ext = '.pdb'
    elif input_file.endswith('.cif'):
        parser = MMCIFParser(QUIET=1)
        ext = '.cif'

    structure = parser.get_structure(pdb_id, input_file)
    chains = structure.get_chains()

    output_files = []
    # Split into chains and save each chain to a new file
    for chain in chains:
        chain_id = chain.get_id()
        outfile = os.path.join(temp_dir, f"{pdb_id}_{chain_id}{ext}")
        io_w_no_h = PDBIO()
        io_w_no_h.set_structure(chain)
        io_w_no_h.save(outfile, ChainSelect(chain_id))
        output_files.append(outfile)

    return output_files

def filter_residues_by_plddt(input_file, plddt_threshold, minimum_plddt):
    """ Filter out residues with sequences of low pLDDT scores """

    temp_dir = tempfile.mkdtemp()
    filename = os.path.basename(input_file)
    output_file = os.path.join(temp_dir, filename)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('structure', input_file)

    residues_to_remove = set()

    for model in structure:
        for chain in model:
            low_plddt_count = 0
            residues_to_check = []
            for residue in chain:
                if 'CA' in residue:
                    plddt = residue['CA'].get_bfactor()
                    print(plddt, plddt_threshold, low_plddt_count)
                    if plddt < plddt_threshold:
                        low_plddt_count += 1
                        residues_to_check.append(residue)
                    else:
                        if low_plddt_count >= minimum_plddt:
                            residues_to_remove.update(residues_to_check)
                        low_plddt_count = 0
                        residues_to_check = []
            if low_plddt_count >= minimum_plddt:
                residues_to_remove.update(residues_to_check)

    io = PDBIO()
    io.set_structure(structure)
    io.save(output_file, ResidueFilterSelect(residues_to_remove))

    return output_file

def preprocess(infile, output_dir, remove_hoh=True, remove_hetatm=True, split_chains=False, alphafold=False, plddt_threshold=65, minimum_plddt=5):

    log_file = os.path.join(output_dir, 'preprocessing_log.txt')
    with open(log_file, 'w') as f:
        f.write(f'Input file: {infile}\n')
        f.write(f'Remove water: {remove_hoh}\n')
        f.write(f'Remove HETATM: {remove_hetatm}\n')
        f.write(f'Split chains: {split_chains}\n')

    preprocessed_pdb = infile

    if remove_hetatm:
        preprocessed_pdb = preprocess_remove_hoh(preprocessed_pdb)
    
    if remove_hetatm:
        preprocessed_pdb = preprocess_remove_hetatm(preprocessed_pdb)

    if split_chains:
        preprocessed_pdb = preprocess_split_chains(preprocessed_pdb)

    if type(preprocessed_pdb) == str:
        preprocessed_pdb = [preprocessed_pdb]

    for file in preprocessed_pdb:
        if alphafold:
            file = filter_residues_by_plddt(file, plddt_threshold, minimum_plddt)
        file_name = os.path.basename(file)
        shutil.copy(file, os.path.join(output_dir, "processed_" + file_name))

    return preprocessed_pdb
