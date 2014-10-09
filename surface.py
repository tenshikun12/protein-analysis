'''
Implements Surface class

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
from math import sin, cos, asin, acos
from point import Point

class Surface(Point):
    '''
    Each of the sections on an atoms surface
    '''
    def __init__(self, polar_coords, atom, size):
        '''
        Surface object constructor
        '''
        Point.__init__(self, None)
        self.radius, self.theta, self.phi = polar_coords
        self.atom = atom
        self.size = size

    _pos_x = None
    @property
    def pos_x(self):
        '''
        Returns cartesian X
        '''
        return self.atom.pos_x+self.radius*sin(self.theta)*cos(self.phi)

    _pos_y = None
    @property
    def pos_y(self):
        '''
        Returns cartesian Y
        '''
        return self.atom.pos_y+self.radius*sin(self.theta)*sin(self.phi)

    _pos_z = None
    @property
    def pos_z(self):
        '''
        Returns cartesian Z
        '''
        return self.atom.pos_z+self.radius*cos(self.theta)

    def __add__(self, other):
        '''
        Sum two surface areas
        '''
        assert isinstance(other, Surface)
        return self.size+other.size

    def __str__(self):
        '''
        Print the properties of the surface
        '''
        return '('+str(self.radius)+','+str(self.theta)+','+\
                str(self.phi)+'): '+str(self.size)+' in '+str(self.atom)

    def __lt__(self, other):
        '''
        Define how to compare two surfaces
        '''
        if abs(self.size-other.size) < 1.E-10:
            if abs(self.radius-other.radius) < 1.E-10:
                if abs(self.theta-other.theta) < 1.E-10:
                    return self.phi < other.phi
                else:
                    return self.theta < other.theta
            else:
                return self.radius < other.radius
        else:
            return self.size < other.size

######### Related to surface but not member of the class

def limit_theta(layer, total_layers):
    '''
    Gives the (upper) limit polar angle for layer k out of a total of m
    Division made for the whole sphere!!
    '''
    assert layer <= total_layers, 'k needs to be \
        less than or equal to m'
    if layer == total_layers:
        return acos(-1.)
    else:
        return acos(1-2.*((2*layer+1.)/(2*total_layers+1))**2)

def center_theta(layer, total_layers):
    '''
    Gives the polar angle to the center of layer layer out of total_layers
    Division made for the whole sphere!!
    '''
    assert layer <= total_layers,  'layer needs to be \
        less than or equal to total_layers'
    if layer == 0:
        return 0.
    else:
        return acos(1-2.*(2.*layer/(2*total_layers+1))**2)

def limit_latitude(layer, total_layers):
    '''
    Gives the (lower) limit latitude angle for layer layer out of total_layers
    Division made for half a sphere!!
    '''
    assert layer <= total_layers, 'layer needs to be \
        less than or equal to total_layers'
    return asin(1-((2*layer+1.)/(2*total_layers+1))**2)

def center_latitude(layer, total_layers):
    '''
    Gives the latitude angle to the center of layer layer out of total_layers.
    Division made for half a sphere!!
    '''
    assert layer <= total_layers,  'layer needs to be \
        less than or equal to total_layers'
    if layer == 0:
        return acos(-1.)/2
    else:
        return asin(1-(2.*layer/(2*total_layers+1))**2)

def center_phi(section, layer):
    '''
    Gives the azimuth angle to the center of section section in layer layer
    '''
    assert section <= 8*layer,  'section section \
        needs to be less than or equal to 8*layer'
    return acos(-1)*(2*section+1.)/(8*layer)

def limit_phi(section, layer):
    '''
    Gives the (upper) limit azimuth angle of section section in layer layer
    '''
    assert section <= 8*layer,  'section section needs to be less \
        than or equal to 8*layer'
    return acos(-1)*2.*section/(8*layer)

def compute_surfaces(atom, total_layers, solvent_radius):
    '''
    Divides the surface of the atom in 1+m*(m+1)/2 sections
    '''
    pi = acos(-1.)
    surf_list = set()
    radius = atom.radius+solvent_radius
    size = 2*pi*(1-cos(atom.limit_theta(0, total_layers)))*radius**2
    surf_list.add(Surface((radius, 0., 0.), atom, size))
    for layer in range(1, total_layers+1):
        theta_increment = cos(atom.limit_theta(layer-1, total_layers))-\
            cos(atom.limit_theta(layer, total_layers))
        for section in range(0, 8*layer):
            phi_increment = atom.limit_phi(section+1, layer)-\
                atom.limit_phi(section, layer)
            if phi_increment < 0:
                phi_increment += 2.*pi
            size = phi_increment*theta_increment*radius**2.
            surf_list.add(Surface((radius, \
                atom.center_theta(layer, total_layers),\
                atom.center_phi(section, layer)), atom, size))
    return surf_list

def compute_surfaces_csu(atom, nko, solvent_radius):
    '''
    Computes the surfaces as done in the csu code
    '''
    pi = acos(-1)
    radius = atom.radius+solvent_radius
    surf_list = []
    fibonacci_a = 0
    fibonacci_b = 1
    for i in range(nko):
        fibonacci_aux = fibonacci_a+fibonacci_b
        fibonacci_a = fibonacci_b
        fibonacci_b = fibonacci_aux
    size = 4*pi*radius/fibonacci_b
    phi_aux = 0
    for k in range(fibonacci_b):
        surf_x = float(k+1)/float(fibonacci_b)
        phi_aux += fibonacci_a
        if phi_aux > fibonacci_b:
            phi_aux -= fibonacci_b
        surf_y = float(phi_aux)/float(fibonacci_b)
        phi = 2*pi*surf_y
        theta = acos(1.-2*surf_x)
        surf_list.append(Surface((radius, \
                theta, phi), atom, size))
    return surf_list

def compute_surfaces_z_simmetry(atom, total_layers, solvent_radius):
    '''
    Divides the surface of the upper half of the sphere in $1+4*(total_layers+1)*total_layers$. The divisions in the lower half are symmetric to this. Returns a set with allsections (both upper and lower halfs).
    '''
    pi = acos(-1.)
    surf_list = set()
    radius = atom.radius+solvent_radius
    size = 2*pi*(1.-sin(limit_latitude(0, total_layers)))*radius**2
    surf_list.add(Surface((radius, 0., 0.), atom, size))
    for layer in range(1, total_layers+1):
        theta_increment = sin(limit_latitude(layer-1, total_layers))-\
            sin(limit_latitude(layer, total_layers))
        for section in range(0, 8*layer):
            phi_increment = limit_phi(section+1, layer)-\
                limit_phi(section, layer)
            if phi_increment < 0:
                phi_increment += 2.*pi
            size = phi_increment*theta_increment*radius**2.
            surf_list.add(Surface((radius, \
                pi/2.-center_latitude(layer, total_layers),\
                center_phi(section, layer)), atom, size))
    return surf_list
