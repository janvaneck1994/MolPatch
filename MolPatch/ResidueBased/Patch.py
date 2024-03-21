class ResiduePatch():
    """
    Parse pisite files to dict

    Attributes
    ----------
    ids : list
        residue ids
    file: str
        pdb file
    dssp: dict
        dssp dict
    dssp_dict_keys: dict
        dssp dict keys
    """

    def __init__(self, ids, dssp, dssp_dict_keys):
        self.ids = ids
        self.dssp = dssp
        self.dssp_dict_keys = dssp_dict_keys

    def get_ids(self):
        return self.ids

    def size(self):
        return sum(self.dssp[i[-2:]][2] for i in self.ids)

    def residues(self):
        return [self.dssp[i[-2:]][0] for i in self.ids]

    def patch_length(self):
        return len(self.ids)

    def residue_on_surface(self):
        """
        Get all residues accessible from the surface

        Return
        ------
        list
            residue ids
        """
        residue_list = []
        for x in self.dssp_dict_keys:
            if x[1][0] != ' ' or  x[1][2] != ' ':
                continue
            if self.dssp[x][2] > 0:
                residue_list.append(x)
        return residue_list
