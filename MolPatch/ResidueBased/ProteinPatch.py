from Bio.PDB import PDBParser, MMCIFParser, Selection
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from ResidueBased.Patch import ResiduePatch
from scipy.spatial import KDTree
from Bio.SeqUtils import seq1
import networkx as nx
import numpy as np

class ProteinPatch():

    def __init__(self, id, file, residues_in_patch, r=1.25, msms='msms -density 1.5'):
        if file.endswith('.pdb'):
            parser = PDBParser()    
        elif file.endswith('.cif'):
            parser = MMCIFParser()
        structure = parser.get_structure(id, file)
        self.model = structure[0]
        self.r = r
        if not msms.startswith('msms'):
            msms = 'msms ' + msms
        self.msms = msms
        self.residues_in_patch = residues_in_patch
        self.dssp_dict, self.dssp_dict_keys =  dssp_dict_from_pdb_file(file)
        self.G = self.dot_cloud_graph()
        self.G = self.patch_network(self.G)
        self.patches = self.create_patches()

    def dot_cloud_graph(self):
        """
        Create a dotted cloud of the proteins surface and label the dots
        hydrophobic or hydrophilic
        ...

        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic or not
        """
        surface_points = get_surface(self.model, MSMS=self.msms)
        residue_list = [r for r in Selection.unfold_entities(self.model, "R") if seq1(r.get_resname()) != 'X']
        center_vectices = [self._sidechain_center(r.get_atoms()) for r in residue_list]

        T = KDTree(center_vectices)
		
        closest_residues = T.query(surface_points, k=1)[1]

        G = nx.Graph()
        for node, coordinates in enumerate(surface_points):
            G.add_node(node)
            G.nodes[node]['selected'] = 0
            closest_residue = residue_list[closest_residues[node]]
            if seq1(closest_residue.get_resname()) in self.residues_in_patch:
                G.nodes[node]['selected'] = 1
            G.nodes[node]['surface_vector_pos'] = coordinates
            G.nodes[node]['closest_residue_id'] = closest_residue.get_full_id()
            G.nodes[node]['closest_residue_aa'] = seq1(closest_residue.get_resname())
        return G

    def patch_network(self, G):
        """
        Create a edges between hydrophobic nodes if they are within r
        hydrophobic or hydrophilic
        ...

        Attributes
        ----------
        G : networkx.Graph
            Graph with nodes labeled as hydrophobic or not

        Return
        ------
        networkx.Graph
            Graph with nodes labeled as hydrophobic and edges between nodes within range r
        """
        node_list = [i for i in G.nodes if G.nodes[i]['selected']]
        x = [G.nodes[i]['surface_vector_pos'] for i in G.nodes if G.nodes[i]['selected']]
        T = KDTree(x)
        pairs = T.query_pairs(self.r)
        G.add_edges_from([(node_list[x[0]],node_list[x[1]]) for x in pairs])
        return G

    def largest_patch(self):
        """
        get largest patch

        Return
        ------
        networkx.Graph
            Graph of the largest patch
        """
        largest_patch = max(self.patches, key=(lambda x: x.size()))
        return largest_patch

    def create_patches(self):
        """
        Get the components of the graph

        Return
        ------
        list
            list of Graph components where an item is a patch
        """
        patched_G = self.patch_network(self.G)
        components = [patched_G.subgraph(c) for c in nx.connected_components(patched_G)] #nx.connected_component_subgraphs(patched_G)
        patch_dict = []
        for component in components:
            if len(component.nodes) <= 1:
                continue
            residue_ids_in_patch = list(set([component.nodes[i]['closest_residue_id'] for i in component.nodes]))
            patch_dict.append(ResiduePatch(residue_ids_in_patch, self.dssp_dict, self.dssp_dict_keys))

        return sorted(patch_dict, key=(lambda x: x.size()), reverse=True)

    def _sidechain_center(self, atoms):
        vectors = [atom.get_vector().get_array() for atom in atoms]
        center = np.array(vectors).mean(axis=0)
        return center
