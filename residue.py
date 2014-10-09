'''
Implements Residue class

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

from point import Point

class Residue(Point):
    '''
    Defines an atom and stores all its properties
    '''
    def __init__(self, num, name):
        '''
        Residue creator
        Receives the residue number in the protein and the name of the residue (3-letter code)
        '''
        Point.__init__(self, None)
        self.num = num
        self.name = name
        self.atoms = set()
        self.structure = None

    @property
    def pos_x(self):
        '''
        Returns X of the (geometrical) center of the residue, taking into account only the heavy atoms.
        '''
        heavy_atoms = [x for x in self.heavy_atoms()]
        if len(heavy_atoms) == 0:
            ret = None
        else:
            pos = 0
            count = 0
            for atom in self.heavy_atoms():
                pos += atom.pos_x
                count += 1
            ret = pos/count
        return ret

    @property
    def pos_y(self):
        '''
        Returns Y of the (geometrical) center of the residue, taking into account only the heavy atoms.
        '''
        heavy_atoms = [x for x in self.heavy_atoms()]
        if len(heavy_atoms) == 0:
            ret = None
        else:
            pos = 0
            count = 0
            for atom in self.heavy_atoms():
                pos += atom.pos_y
                count += 1
            ret = pos/count
        return ret

    _pos_z = None
    @property
    def pos_z(self):
        '''
        Returns Z of the (geometrical) center of the residue, taking into account only the heavy atoms.
        '''
        heavy_atoms = [x for x in self.heavy_atoms()]
        if len(heavy_atoms) == 0:
            ret = None
        else:
            pos = 0
            count = 0
            for atom in self.heavy_atoms():
                pos += atom.pos_z
                count += 1
            ret = pos/count
        return ret
    
    _radius = None
    @property
    def radius(self):
        '''
        Returns maximum dimension in the residue.
        '''
        if self._radius == None:
            min_x = 1000
            max_x = -1000
            min_y = 1000
            max_y = -1000
            min_z = 1000
            max_z = -1000
            for atom in self.heavy_atoms():
                if atom.pos_x > max_x:
                    max_x = atom.pos_x
                if atom.pos_y > max_y:
                    max_y = atom.pos_y
                if atom.pos_z > max_z:
                    max_z = atom.pos_z
                if atom.pos_x < min_x:
                    min_x = atom.pos_x
                if atom.pos_y < min_y:
                    min_y = atom.pos_y
                if atom.pos_z < min_z:
                    min_z = atom.pos_z
            self._radius = max([max_x-min_x, \
                    max_y-min_y, max_z-min_z])
        return self._radius
    
    def heavy_atoms(self):
        '''
        Yields the heavy atoms in the residue
        '''
        for atom in self.atoms:
            if atom.heavy:
                yield atom

    def get_atom(self, name):
        '''
        Returns the atom in the residue with the given name
        '''
        lst = [x for x in self.atoms if x.name == name]
        assert len(lst) < 2, "Apparently there is more than one atom \
                with name '%r'" % name
        if len(lst) == 1:
            return lst[0]
        elif len(lst) == 0:
            return None
        
    def add_atom(self, atom):
        '''
        Adds the given atom to the residue
        Throws an exception if there is already an atom in the residue with that name
        '''
        assert self.get_atom(atom.name) == None, "Residue already \
                contains a '%r'" % atom.name
        self.atoms.add(atom)
        atom.set_residue(self)
