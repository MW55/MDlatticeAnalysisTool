from Bio.PDB.Atom import Atom
from Bio.PDB.PDBParser import PDBParser
import numpy as np
import itertools

class Enviroment(object):
    def __init__(self, protein_pdb_path, polymer_full_poses_path,
            meshsize = np.array([3, 3, 3])):
        self.prot_path = protein_pdb_path
        self.poly_path = polymer_full_poses_path
        self._meshsize = meshsize
        self._protein_struc = self._read_protein(self.prot_path)
        self._poly_poses = self._read_polymer_poses(self.poly_path)
        self.protein_coords = self._get_protein_coords(self._protein_struc)
        self.poly_poses_coords = self._get_poly_coords(self._poly_poses)
        self.min_point = self.min_point(self.poly_poses_coords,
                self.protein_coords)
        self.max_point = self.max_point(self.poly_poses_coords,
                self.protein_coords)
        self._residues = self._filter_residues(self._protein_struc)
        self.residue_list = self._create_residue_list(self._residues)
        self.res_id_dict = self._create_residue_id_dict(self._residues)
        self.geometric_center = self._get_geometric_center(self.protein_coords)

    def _read_protein(self, prot_path):
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        return(parser.get_structure(prot_path, prot_path))

    def _read_polymer_poses(self, poly_path):
        parser = PDBParser(PERMISSIVE=1, QUIET=True)
        structure = parser.get_structure(poly_path, poly_path)
        poses = [model for model in structure]
        return poses

    def _get_protein_coords(self, protein_struc):
        protein_coords = [atom.get_coord() for atom in protein_struc.get_atoms()]
        return protein_coords

    def _get_poly_coords(self, poly_poses):
        poses_coords = []
        for model in poly_poses:
            model_coords = []
            for atom in model.get_atoms():
                model_coords.append(atom.get_coord())
            poses_coords.append(model_coords)
        return poses_coords

    def _get_geometric_center(self, protein_coords):
        return(sum(protein_coords)/len(protein_coords))

    # Calculate the lowest cartesian coordinates (rounded down) of all poses.
    def min_point(self, poly_poses_coords, protein_coords):
        minx = min([protein_coords[atom][0] for atom in
            range(len(protein_coords))])
        miny = min([protein_coords[atom][1] for atom in
            range(len(protein_coords))])
        minz = min([protein_coords[atom][2] for atom in
            range(len(protein_coords))])
        for i in range(len(poly_poses_coords)):
            for j in range(len(poly_poses_coords[i])):
                if poly_poses_coords[i][j][0] < minx:
                    minx = poly_poses_coords[i][j][0]
                if poly_poses_coords[i][j][1] < miny:
                    miny = poly_poses_coords[i][j][1]
                if poly_poses_coords[i][j][2] < minz:
                    minz = poly_poses_coords[i][j][2]
        return(np.array([np.floor(minx), np.floor(miny), np.floor(minz)]))

    # Calculate the highest cartesian coordinates (rounded up) of all poses.
    def max_point(self, poly_poses_coords, protein_coords):
        maxx = max([protein_coords[atom][0] for atom in
            range(len(protein_coords))])
        maxy = max([protein_coords[atom][1] for atom in
            range(len(protein_coords))])
        maxz = max([protein_coords[atom][2] for atom in
            range(len(protein_coords))])
        for i in range(len(poly_poses_coords)):
            for j in range(len(poly_poses_coords[i])):
                if poly_poses_coords[i][j][0] > maxx:
                    maxx = poly_poses_coords[i][j][0]
                if poly_poses_coords[i][j][1] > maxy:
                    maxy = poly_poses_coords[i][j][1]
                if poly_poses_coords[i][j][2] > maxz:
                    maxz = poly_poses_coords[i][j][2]
        return(np.array([np.ceil(maxx), np.ceil(maxy), np.ceil(maxz)]))

    def _filter_residues(self, protein_struc):
        filter_hetatm = lambda x: x.id[0]==' '
        return(filter(filter_hetatm, protein_struc.get_residues()))

    # Creates a list of the residues with the coordinates of all atoms
    # and the corresponding residue id as values.
    def _create_residue_list(self, residues):
        resi_coord_id_list = []
        for residue in residues:
            resi_coord_id_list.append([np.append(atom.get_coord(), residue.id[1])
                for atom in residue.get_list()])
        return list(itertools.chain.from_iterable(resi_coord_id_list))

    # Creates a dict containing the residue ids and the corresponding
    # residues.
    def _create_residue_id_dict(self, residues):
        residdict = {}
        for residue in residues:
            residdict[residue.id[1]] = residue.resname
        return residdict

# The box generation is not in use at the moment, instead the box created
# by epitopsy is imported, to allow for comparisons between epitopsy
# results and mdl results.

    # def _create_box(self):
    #     box = Box(self)
    #     return box

# class Box(object):
#     def __init__(self, enviroment):
#         self.offset = self.offset(enviroment.min_point)
#         self.meshsize = enviroment._meshsize
#         self.box_dim = self._calc_box_dim(enviroment.max_point,
#                enviroment.min_point, self.meshsize)
#         self.meshes_amount = self._number_of_meshes(self.meshsize,
#                   self.box_dim)

#     # is used to translate the atomic coordinates to the coordinates
      # of the box and vice versa.
#     def offset(self, min_point):
#         return abs(min_point)

#     def _calc_box_dim(self, max_point, min_point, meshsize):
#         box_dim = np.zeros(3)
#         box_dim[0] = abs(max_point[0] - min_point[0])
#         box_dim[1] = abs(max_point[1] - min_point[1])
#         box_dim[2] = abs(max_point[2] - min_point[2])
#         return np.floor(box_dim / self.meshsize)

#     def _number_of_meshes(self, meshsize, box_dim):
#         return np.prod(box_dim)

#     def box_to_real_space(self, grid_coord):
#         return grid_coord * self.meshsize - self.offset

#     def real_to_box_space(self, atom_coord):
#         return np.around((atom_coord + self.offset) / self.meshsize)
