'''
Implements Atom class

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
import surface

class Atom(Point):
    '''
    Defines an atom and stores all its properties
    '''
    def __init__(self, num, name, position):
        '''
        Atom creator
        Receives the index in the protein,  the name (CA, OE2...) and the position coordinates
        '''
        Point.__init__(self, position)
        self.num = num
        self.name = name
        self.residue = None
        self.heavy = True
        self._surfaces = {}
        if name[0] in ['H', '1', '2']:
            self.heavy = False

    def mass(self, refactor = 1.66E-27):
        '''
        Returns the mass of the atom.
        If refactor is 1,  mass is returned in atomic mass units.
        The default returns it in kg.
        '''
        ret = None
        elem = self.name[0]
        if elem == 'C':
            ret = 12
        elif elem == 'O':
            ret = 16
        elif elem == 'N':
            ret = 14
        elif elem == 'S':
            ret = 32
        return ret*refactor

    def __str__(self):
        '''
        Returns a string with the characteristics of the atom
        '''
        return '%(resName)s%(resId)d:%(name)s' % {'resName':self.residue.name, \
                'resId':self.residue.num, 'name':self.name}

    def __lt__(self, other):
        '''
        Custom comparison:
        * a1 is smaller than a2 if the residue of a1 has a lower index.
        * if in the same residue, the index of the atoms is compared.
        '''
        if self.residue.num == other.residue.num:
            return self.num < other.num
        else:
            return self.residue.num < other.residue.num

    def set_residue(self, residue):
        '''
        Residue setter: Assigns a residue to the atom
        '''
        self.residue = residue

    def get_surfaces(self, total_layers, solvent_radius, \
            compute_function = surface.compute_surfaces):
        '''
        Returns the surfaces of the atom.
        Exists because once computed you don't need to compute them again
        '''
        try:
            ret = self._surfaces[(total_layers, solvent_radius)]
        except KeyError:
            ret = compute_function(self, total_layers, solvent_radius)
            self._surfaces[(total_layers, solvent_radius)] = ret
        return ret

    _csu_radius = {'N': 1.7, 'O': 1.5, 'C': 1.9, 'S': 1.9}
    _radius = {}
    #Backbone
    _radius['N'] = 1.64
    _radius['CA'] = 1.88
    _radius['C'] = 1.61
    _radius['O'] = 1.42
    _radius['OC1'] = 1.42
    _radius['OXT'] = 1.46
    _radius['OC2'] = 1.46
    #Beta
    _radius['CB'] = {'ALA': 1.88, 'VAL': 1.88, 'PRO': 1.88, 'LEU': 1.88, \
            'ILE': 1.88, 'PHE': 1.88, 'MET': 1.88, 'ASP': 1.88, 'GLU': 1.88, \
            'LYS': 1.88, 'ARG': 1.88, 'SER': 1.88, 'THR': 1.88, 'TYR': 1.88, \
            'HIS': 1.88, 'CYS': 1.88, 'ASN': 1.88, 'GLN': 1.88, 'TRP': 1.88}
    #Gamma
    _radius['CG'] = {'PRO': 1.88, 'LEU': 1.88, 'PHE': 1.88, 'MET': 1.88, \
            'ASP': 1.61, 'GLU': 1.88, 'LYS': 1.88, 'ARG': 1.88, 'TYR': 1.88, \
            'HIS': 1.61, 'ASN': 1.61, 'GLN': 1.88, 'TRP': 1.61}
    _radius['CG1'] = {'VAL': 1.88, 'ILE': 1.88}
    _radius['CG2'] = {'VAL': 1.88, 'ILE': 1.88, 'THR': 1.88}
    _radius['OG'] = {'SER': 1.46}
    _radius['OG1'] = {'THR': 1.46}
    _radius['SG'] = {'CYS': 1.77}
    #Delta
    _radius['CD'] = {'GLU': 1.61, 'LYS': 1.88, 'ARG': 1.88, 'GLN': 1.61, \
            'PRO': 1.88}
    _radius['CD1'] = {'LEU': 1.88, 'ILE': 1.88, 'PHE': 1.61, 'TYR': 1.76, \
            'TRP': 1.76}
    _radius['CD2'] = {'LEU': 1.88, 'PHE': 1.76, 'TYR': 1.76, 'HIS': 1.76, \
            'TRP': 1.61}
    _radius['OD1'] = {'ASP': 1.46, 'ASN': 1.42}
    _radius['OD2'] = {'ASP': 1.42}
    _radius['SD'] = {'MET': 1.77}
    _radius['ND1'] = {'HIS': 1.64}
    _radius['ND2'] = {'ASN': 1.64}
    #Epsilon
    _radius['CE'] = {'MET': 1.88, 'LYS': 1.88}
    _radius['CE1'] = {'PHE': 1.76, 'TYR': 1.76, 'HIS': 1.76}
    _radius['CE2'] = {'PHE': 1.76, 'TYR': 1.76, 'TRP': 1.61}
    _radius['CE3'] = {'TRP': 1.76}
    _radius['OE1'] = {'GLU': 1.46, 'GLN': 1.42}
    _radius['OE2'] = {'GLU': 1.42}
    _radius['NE'] = {'ARG': 1.64}
    _radius['NE1'] = {'TRP': 1.64}
    _radius['NE2'] = {'HIS': 1.64, 'GLN': 1.64}
    #Zeta
    _radius['CZ'] = {'PHE': 1.76, 'ARG': 1.61, 'TYR': 1.61}
    _radius['CZ2'] = {'TRP': 1.76}
    _radius['CZ3'] = {'TRP': 1.76}
    _radius['NZ'] = {'LYS': 1.64}
    #Eta (H)
    _radius['CH2'] = {'TRP': 1.76}
    _radius['OH'] = {'TYR': 1.46}
    _radius['NH1'] = {'ARG': 1.64}
    _radius['NH2'] = {'ARG': 1.64}

    _csu_class = {}
    #Backbone
    _csu_class['N'] = {'ALA': 3, 'ARG': 3, 'ASN': 3, 'ASP': 3, 'CYS': 3, \
            'GLU': 3, 'GLN': 3, 'GLY': 3, 'HIS': 3, 'ILE': 3, 'LEU': 3, \
            'LYS': 3, 'MET': 3, 'PHE': 3, 'PRO': 6, 'SER': 3, 'THR': 3, \
            'TRP': 3, 'TYR': 3, 'VAL': 3}
    _csu_class['CA'] = {'ALA': 7, 'ARG': 7, 'ASN': 7, 'ASP': 7, 'CYS': 7, \
            'GLU': 7, 'GLN': 7, 'GLY': 6, 'HIS': 7, 'ILE': 7, 'LEU': 7, \
            'LYS': 7, 'MET': 7, 'PHE': 7, 'PRO': 4, 'SER': 7, 'THR': 7, \
            'TRP': 7, 'TYR': 7, 'VAL': 7}
    _csu_class['C'] = {'ALA': 6, 'ARG': 6, 'ASN': 6, 'ASP': 6, 'CYS': 6, \
            'GLU': 6, 'GLN': 6, 'GLY': 6, 'HIS': 6, 'ILE': 6, 'LEU': 6, \
            'LYS': 6, 'MET': 6, 'PHE': 6, 'PRO': 6, 'SER': 6, 'THR': 6, \
            'TRP': 6, 'TYR': 6, 'VAL': 6}
    _csu_class['O'] = {'ALA': 2, 'ARG': 2, 'ASN': 2, 'ASP': 2, 'CYS': 2, \
            'GLU': 2, 'GLN': 2, 'GLY': 2, 'HIS': 2, 'ILE': 2, 'LEU': 2, \
            'LYS': 2, 'MET': 2, 'PHE': 2, 'PRO': 2, 'SER': 2, 'THR': 2, \
            'TRP': 2, 'TYR': 2, 'VAL': 2}
    _csu_class['OXT'] = {'ALA': 1, 'ARG': 1, 'ASN': 1, 'ASP': 1, 'CYS': 1, \
            'GLU': 1, 'GLN': 1, 'GLY': 1, 'HIS': 1, 'ILE': 1, 'LEU': 1, \
            'LYS': 1, 'MET': 1, 'PHE': 1, 'PRO': 1, 'SER': 1, 'THR': 1, \
            'TRP': 1, 'TYR': 1, 'VAL': 1}
    #Beta
    _csu_class['CB'] = {'ALA': 4, 'ARG': 4, 'ASN': 4, 'ASP': 4, 'CYS': 4, \
            'GLU': 4, 'GLN': 4, 'HIS': 4, 'ILE': 4, 'LEU': 4, 'LYS': 4, \
            'MET': 4, 'PHE': 4, 'PRO': 4, 'SER': 6, 'THR': 6, 'TRP': 4, \
            'TYR': 4, 'VAL': 4}
    #Gamma
    _csu_class['CG'] = {'ARG': 4, 'ASN': 6, 'ASP': 6, 'GLU': 4, 'GLN': 4, \
            'HIS': 5, 'LEU': 4, 'LYS': 4, 'MET': 4, 'PHE': 5, 'PRO': 4, \
            'TRP': 5, 'TYR': 5}
    _csu_class['CG1'] = {'ILE': 4, 'VAL': 4}
    _csu_class['CG2'] = {'ILE': 4, 'THR': 4, 'VAL': 4}
    _csu_class['OG'] = {'SER': 1}
    _csu_class['OG1'] = {'THR': 1}
    _csu_class['SG'] = {'CYS': 6}
    #Delta
    _csu_class['CD'] = {'ARG': 7, 'GLU': 6, 'GLN': 6, 'LYS': 4, 'PRO': 4}
    _csu_class['CD1'] = {'ILE': 4, 'LEU': 4, 'PHE': 5, 'TRP': 5, 'TYR': 5}
    _csu_class['CD2'] = {'HIS': 5, 'LEU': 4, 'PHE': 5, 'TRP': 5, 'TYR': 5}
    _csu_class['ND1'] = {'HIS': 1}
    _csu_class['ND2'] = {'ASN': 3}
    _csu_class['OD1'] = {'ASN': 2, 'ASP': 2}
    _csu_class['OD2'] = {'ASP': 2}
    _csu_class['SD'] = {'MET': 8}
    #Epsilon
    _csu_class['CE'] = {'LYS': 7, 'MET': 4}
    _csu_class['CE1'] = {'HIS': 5, 'PHE': 5, 'TYR': 5}
    _csu_class['CE2'] = {'PHE': 5, 'TRP': 5, 'TYR': 5}
    _csu_class['CE3'] = {'TRP': 5}
    _csu_class['NE'] = {'ARG': 3}
    _csu_class['NE1'] = {'TRP': 3}
    _csu_class['NE2'] = {'GLN': 3, 'HIS': 1}
    _csu_class['OE1'] = {'GLU': 2, 'GLN': 2}
    _csu_class['OE2'] = {'GLU': 2}
    #Zeta
    _csu_class['CZ'] = {'ARG': 6, 'PHE': 5, 'TYR': 5}
    _csu_class['CZ2'] = {'TRP': 5}
    _csu_class['CZ3'] = {'TRP': 5}
    _csu_class['NZ'] = {'LYS': 3}
    #Eta (H)
    _csu_class['CH2'] = {'TRP': 5}
    _csu_class['NH1'] = {'ARG': 3}
    _csu_class['NH2'] = {'ARG': 3}
    _csu_class['OH'] = {'TYR': 1}

    @property
    def csu_class(self):
        '''
        csu_class getter
        '''
        return self._csu_class[self.name][self.residue.name]

    @property
    def radius(self):
        '''
        Returns the Van der Waals radius of each atom
        '''
        #return self._csu_radius[self.name[0]]
        if self.name in ['N', 'CA', 'C', 'O', 'OC1', 'OC2', 'OXT']:
            return self._radius[self.name]
        else:
            return self._radius[self.name][self.residue.name]
