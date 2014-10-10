'''
Implements Protein class

This file is part of protein-analysis

protein-analysis is ree software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

protein-analysis is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with protein-analysis.  If not, see <http://www.gnu.org/licenses/>.
'''

from atom import Atom
from residue import Residue

AMINO_ACID_NAMES = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', \
        'GLN', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', \
        'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']

class Protein(list):
    '''
    Changes list class so that residues can be selected as in most other programs (start with 1, end with len(protein)
    '''
    def __getitem__(self, res_id):
        '''
        Return the residue with index res_id (starting at 1)
        '''
        try:
            return super(Protein, self).__getitem__(res_id-1)
        except IndexError:
            return None

    def __getslice__(self, res_id1, res_id2):
        '''
        Returns all residues between indexes res_id1 (included) and resId2 (not included), starting at 1
        '''
        res_id1 -= 1
        res_id2 -= 1
        return super(Protein, self).__getslice__(res_id1, res_id2)

    def __setitem__(self, index, residue):
        '''
        Sets the value of the index'th element (starting at 1) of the protein to residue
        '''
        assert residue == None or index == residue.num, \
                "Residue has a different number than the index"
        super(Protein, self).__setitem__(index-1, residue)

    def __setslice__(self, index1, index2, residue):
        '''
        Sets the value of all residues from index1 (included) to index2 (not included) to residue, starting at 1.
        Residue must be None (it is not default for security reasons), other values raise exceptions.
        '''
        assert residue == None, "Slice-setting \
                can only be used for initializing to None"
        super(Protein, self).__setslice__(index1-1, index2-1, None)

    def __delitem__(self, index):
        '''
        If index corresponds to the last element of the protein, it deletes the residue. Otherwise, it makes it None (so as to not break the order).
        '''
        if len(self) == index:
            super(Protein, self).__delitem__(index-1)
        else:
            self[index] = None

    def __delslice__(self, index1, index2):
        '''
        If index2 corresponds with the last element, remove residues from index1 until the end of the protein. Otherwise, puts residues between index1 and index2 to None.
        '''
        if len(self) == index2:
            super(Protein, self).__delslice__(index1-1, index2-1)
        else:
            self[index1, index2] = None

    def get_residue(self, res_name):
        '''
        Returns a list of the residues in the protein with name res_name (3-letter code)
        '''
        return [x for x in self if x.name == res_name]

    def append(self, residue):
        '''
        Adds a residue at the end of the protein
        '''
        self.add_residue(residue)

    def insert(self, index, residue):
        '''
        Inserts a residue at position index.
        Can only be used if with the same index as residue.num, otherwise it throws an exception.
        '''
        assert index == residue.num, \
                "Residue has a different number than the index"
        self.add_residue(residue)
            
    def index(self, value):
        '''
        Returns the index of value (starting at 1)
        '''
        return super(Protein, self).index(value)+1

    def add_residue(self, residue):
        '''
        Adds a residue to the protein.
        If another residue with the same number exists, it throws an exception.
        '''
        assert self[residue.num] == None, "Protein already \
                contains a residue number '%d'" % residue.num
        protein_length = len(self)
        while protein_length < residue.num-1:
            super(Protein, self).append(None)
            protein_length = len(self)
        super(Protein, self).insert(residue.num-1, residue)

    def add_atom(self, atom, res_id):
        '''
        Adds atom to the residue corresponding to res_id.
        The residue must exist, and all specifications in Residue.addAtom apply.
        '''
        res = self[res_id]
        assert res != None, "Residue %d does not exist!" % res_id
        res.add_atom(atom)

    def heavy_atoms(self):
        '''
        Yields all heavy atoms in the protein.
        '''
        for res in self:
            if res != None:
                for atom in res.atoms:
                    if atom.heavy:
                        yield atom

    def atoms(self):
        '''
        Yields all atoms in the protein
        '''
        for res in self:
            for atom in res.atoms:
                yield atom

    def set_structure(self, structure_file):
        '''
        Reads secondary structure of the amino acids from an xpm file and associates the structure to each residue.
        '''
        extension = structure_file.split('.')[-1]
        assert extension  ==  'xpm'
        with open(structure_file, 'r') as read_file:
            data = read_file.read().split('\n')
        for i in range(len(self)):
            struc = data[-2-i][1]
            if struc == '~':
                struc = 'C'
            self[i+1].structure = struc

    def get_structure(self):
        '''
        Prints the structure of the protein as a string
        '''
        ret = ''
        for res in self:
            ret += res.structure
        return ret

    def structural_distance(self, other):
        '''
        Computes the distance between two structures using TMalign (TODO)
        '''
        assert isinstance(other, Protein)
        str1 = self.get_structure()
        str2 = other.get_structure()
        ret = abs(len(str1)-len(str2))
        for i in min(len(str1), len(str2)):
            if str1[i] == str2[i]:
                ret += 1
        return ret

    def __init__(self, file_name):
        '''
        Protein class constructor
        Receives the file from which it is to be read. It needs to be in pdb or gro formats.
        Proteins must be ungapped (all resiude numbers must be consecutive).
        Only standard residues are considered.
        '''
        super(Protein, self).__init__()
        ext = file_name.split('.')[-1].lower()
        if ext == 'pdb':
            with open(fileName,'r') as f:
                fileContent=f.read().split('\n')
            for line in fileContent:
                if line[0:6] in ['TER   ','END   ']:
                    return None
                elif line[0:6]=='ATOM  ':
                    index = int(line[6:11])
                    name = line[12:16].strip()
                    resName = line[17:20].strip()
                    if res_name in ['HIP', 'HIE', 'HID']:
                        res_name = 'HIS'
                    assert res_name in AMINO_ACID_NAMES, \
                            "Only standard residues are considered."
                    resId = int(line[22:26])
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if resId > 0:
                        try:
                            self.addResidue(Residue(resId,resName))
                        except AssertionError:
                            pass
                        self.addAtom(Atom(index,name,x,y,z),resId)
            q = [x for x in self if x != None]
            assert None not in self[self.index(q[0]):self.index(q[-1])+1], \
                    "Protein contains gaps."
        elif ext == 'gro':
            with open(file_name, 'r') as read_file:
                file_content = read_file.read().split('\n')
            number_of_atoms = int(file_content[1])
            for index in range(number_of_atoms):
                line = file_content[2+index].split()
                index = int(line[2])
                name = line[1]
                res_name = line[0][-3:]
                res_id = int(line[0][:-3])
                pos_x = float(line[3])*10
                pos_y = float(line[4])*10
                pos_z = float(line[5])*10
                try:
                    self.add_residue(Residue(res_id, res_name))
                except AssertionError:
                    pass
                self.add_atom(Atom(index, name, (pos_x, pos_y, pos_z)), res_id)
        else:
            raise NameError('Unknown file extension')
        return None
