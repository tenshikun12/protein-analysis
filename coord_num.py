'''
Implements functions to generate contact maps and coordination numbers.

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
from math import acos
from protein import Protein
#from residue import Residue
#from atom import Atom
import surface

BOND_TYPE = {\
    1: {1: 'HB', 2: 'HB', 3: 'HB', 4: 'HH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}, \
    2: {1: 'HB', 2: 'AA', 3: 'HB', 4: 'HH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}, \
    3: {1: 'HB', 2: 'HB', 3: 'DD', 4: 'HH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}, \
    4: {1: 'HH', 2: 'HH', 3: 'HH', 4: 'PH', 5: 'PH', 6: 'PH', 7: 'PH', \
        8: 'PH'}, \
    5: {1: 'O', 2: 'O', 3: 'O', 4: 'PH', 5: 'AR', 6: 'O', 7: 'O', 8: 'O'}, \
    6: {1: 'O', 2: 'O', 3: 'O', 4: 'PH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}, \
    7: {1: 'O', 2: 'O', 3: 'O', 4: 'PH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}, \
    8: {1: 'O', 2: 'O', 3: 'O', 4: 'PH', 5: 'O', 6: 'O', 7: 'O', 8: 'O'}}

def avg_coord_num(res2res):
    '''
    Input: dictionary containing all residues and a set of the residues in contact with it (excluding itself, including repetitions)
    Output: Average coordination number: Sum of the number of contacts of each residue normalized to the chain length
    '''
#    print 'DEBUG: avgCoordNum start'
    coord_num = {}
    for k in res2res.keys():
        coord_num[k] = len(res2res[k])
    ret = sum(coord_num.values())/float(len(coord_num.keys()))
#    print 'DEBUG: avgCoordNum end. Returning '+str(ret)
    return ret

def marek_avg_coord_num(protein, alpha = (26/7.)**(1/6.)):
    '''
    Default value is from Marek's original code
    Input: protein as read froma a PDB file
    Output: A tuple containing the number of contacts (j-i>2) and
		    the average coordination number (see avgCoordNum)
    '''
#    print 'DEBUG: marekAvgCoordNum start'
    res2res = {}
    couples = set()
    kcont = 0
    for res1 in [x for x in protein if x != None]:
        res2res[res1] = set()
        for res2 in [x for x in protein if x != res1 and x != None]:
            icc = 0
            for atom1 in res1.heavy_atoms():
                for atom2 in res2.heavy_atoms():
                    if atom1.dist(atom2)  <=  \
                        (atom1.radius+atom2.radius)*alpha:
                        icc += 1
            if icc > 0:
                res2res[res1].add(res2)
                if protein.index(res2)-protein.index(res1)>2:
                    kcont += 1
                couples.add((res1.num, res2.num))
    ret = kcont, avg_coord_num(res2res), res2res
#    print 'DEBUG: marekAvgCoordNum end. Returning '+str(ret)
    return ret

def csu_contact_map(protein, solvent_radius = 1.4, total_layers = 5, \
        compute_function = None):
    '''
    Surface division of spheres done with equal area (Sobolev et al 1996 Proteins 25:120-9)
    Default values are radius of water molecule and best fit to original CSU
    Input: protein as read froma a PDB file
    Output: A tuple containing the number of contacts (j-i>2) and 
		    the average coordination number (see avgCoordNum)
    '''
#    print 'DEBUG: csuAvgCoordNum start'
    atom2atom = {}
    heavy_atoms = [x for x in protein.heavy_atoms()]
    for atom1 in heavy_atoms:
        near = []
        for atom2 in [x for x in heavy_atoms if x.residue != atom1.residue]:
        #for atom2 in [x for x in heavy_atoms if x != atom1]:
            if atom1.dist(atom2) < atom1.radius+atom2.radius+2*solvent_radius:
                near.append(atom2)
        surfaces1 = atom1.get_surfaces(total_layers, solvent_radius, \
                compute_function)
        for surf in surfaces1:
            if len(near)>0:
                atom2 = min(near, key=(lambda x: surf.dist(x)))
                if surf.dist(atom2) < solvent_radius+atom2.radius:
                    try:
                        atom2atom[atom2].add(surf.atom)
                    except KeyError:
                        atom2atom[atom2] = set([surf.atom])
    ret = atom2atom
#    print 'DEBUG: csuAvgCoordNum end. Returning '+str(ret)
    return ret

def really_csu_contact_map(protein, solvent_radius = 1.4, total_layers = 14, \
        compute_function = None):
    '''
    Surface division of spheres done with equal area (Sobolev et al 1996 Proteins 25:120-9)
    Default values are radius of water molecule and best fit to original CSU
    Input: protein as read froma a PDB file
    Output: A tuple containing the number of contacts (j-i>2) and 
		    the average coordination number (see avgCoordNum)
    '''
#    print 'DEBUG: csuAvgCoordNum start'
    atom2atom = {}
    accessible_area = {}
    heavy_atoms = [x for x in protein.heavy_atoms()]
    for atom1 in heavy_atoms:
        atom_to_surface = {}
        atom2atom[atom1] = {}
        sphere_radius = atom1.radius + solvent_radius
        contacting = []
        for residue in [x for x in protein if x != None]:
            if atom1.residue == residue or \
                residue.dist(atom1) <= sphere_radius + solvent_radius \
                + residue.radius:
                for atom2 in residue.heavy_atoms():
                    if atom1 != atom2 and \
                        atom2.dist(atom1) <= sphere_radius + \
                        solvent_radius + residue.radius:
                        contacting.append(atom2)
        #number_contacts = len(contacting)
        contacting.sort(key = lambda x: atom1.dist(x))
        surfaces1 = atom1.get_surfaces(total_layers, solvent_radius, \
                surface.compute_surfaces_csu)
        section_size = 4*acos(-1)*sphere_radius*sphere_radius/ \
                float(len(surfaces1))
        for surf in surfaces1:
            cont = True
            for atom2 in contacting:
                if cont and surf.dist(atom2) < atom2.radius + solvent_radius:
                    cont = False
                    try:
                        atom_to_surface[atom2].add(surf)
                    except KeyError:
                        atom_to_surface[atom2] = set([surf])
        contacting_sections = sum(len(x) for x in atom_to_surface.values())
        accessible_area[atom1] = section_size* \
                (len(surfaces1)-contacting_sections)
        #ikr(iq) = len(atom_to_surface[atom2])
        for atom2 in atom_to_surface.keys():
            if atom2.residue != atom1.residue:
                atom2atom[atom1][atom2] = len(atom_to_surface[atom2])* \
                        section_size
    ret = atom2atom
#    print 'DEBUG: csuAvgCoordNum end. Returning '+str(ret)
    return ret

def csu_avg_coord_num(protein, solvent_radius = 1.4, total_layers = 5, \
        compute_function = surface.compute_surfaces_csu):
    '''
    Surface division of spheres done with equal area (Sobolev et al 1996 Proteins 25:120-9)
    Default values are radius of water molecule and best fit to original CSU
    Input: protein as read froma a PDB file
    Output: A tuple containing the number of contacts (j-i>2) and 
		    the average coordination number (see avgCoordNum)
    '''
#    print 'DEBUG: csuAvgCoordNum start'
    atom2atom = really_csu_contact_map(protein, solvent_radius, total_layers, \
            compute_function)
    res2res = {}
    for atom1 in atom2atom.keys():
        res1 = atom1.residue
        for atom2 in atom2atom[atom1]:
            res2 = atom2.residue
            if res2 != res1:
                try:
                    res2res[res1].add(res2)
                except KeyError:
                    res2res[res1] = set([res2])
    ret = sum([len(res2res[x]) for x in res2res.keys()])/2, \
            avg_coord_num(res2res)
#    print 'DEBUG: csuAvgCoordNum end. Returning '+str(ret)
    return ret

#def picContactNumber(protein):
#    raise NotImplementedError

if __name__  == "__main__":
    import sys
    if len(sys.argv) == 1 or '-h' in sys.argv:
        print "Usage: python coordNum.py pdbFileName1 [pdbFileName2 ...]"
    else:
        #print "Protein name\tContact Number (Marek)\tAvg Coord Number \
        #    (Marek)\tContact Number (CSU)\tAvg Coord Number (CSU)"
        for f in sys.argv[1:]:
            p = Protein(f)
            name = f.split('/')[-2].split('.')[0]
            marek = marek_avg_coord_num(p)
            csu = csu_avg_coord_num(p, total_layers = 14)
            print '%(n)5s %(ncm)3d %(cnm)6.3f %(ncc)3d %(cnc)6.3f' % \
                    {'n': name, 'ncm': marek[0], 'cnm': marek[1], \
                    'ncc': csu[0], 'cnc': csu[1]}
