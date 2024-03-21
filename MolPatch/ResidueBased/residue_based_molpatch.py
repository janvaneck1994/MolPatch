import pandas as pd
import os
from ResidueBased.ProteinPatch import ProteinPatch

def calculate_patches(infile, residues, plot, output_dir):
    """ Calculate patches on a protein surface """

    filename = os.path.basename(infile)
    outfile = filename.split('.')[0] + '.csv'
    pdb_id = filename.split('.')[0].split('_')[0]

    proteinPatches = ProteinPatch(pdb_id,infile,residues)
    patches = proteinPatches.patches
    result_dict = {'patch_rank':[], 'protein_id':[], 'residue_id':[], 'chain':[], 'residue_type':[], 'patch_size':[]}
    for i,patch in enumerate(patches):
        residues_in_patch = patch.get_ids()
        residues = patch.residues()
        
        for j, residue_in_patch in enumerate(residues_in_patch):
            result_dict['residue_id'].append(residue_in_patch[-2:][1][1])
            result_dict['chain'].append(residue_in_patch[-2:][0])
            result_dict['patch_size'].append(patch.size())
            result_dict['patch_rank'].append(i)
            result_dict['residue_type'].append(residues[j])
            result_dict['protein_id'].append(pdb_id)

    pd.DataFrame(result_dict).to_csv(os.path.join(output_dir, "patches_" + outfile), index=False)

    if plot:
        outfig = infile.split('.')[0] + '.png'
        proteinPatches.plot_largest_patches(outfig)
    