'''
Implements Point class

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
from math import sqrt

class Point(object):
    '''
    All objects with a position get the same distance calculation
    '''
    def __init__(self, position):
        '''
        Point constructor
        '''
        if position != None:
            self.pos_x, self.pos_y, self.pos_z = position

    _pos_x = None
    @property
    def pos_x(self):
        '''
        Returns X
        '''
        return self._pos_x
    @pos_x.setter
    def pos_x(self, value):
        '''
        Sets X
        '''
        self._pos_x = value

    _pos_y = None
    @property
    def pos_y(self):
        '''
        Returns Y
        '''
        return self._pos_y
    @pos_y.setter
    def pos_y(self, value):
        '''
        Sets Y
        '''
        self._pos_y = value

    _pos_z = None
    @property
    def pos_z(self):
        '''
        Returns Z
        '''
        return self._pos_z
    @pos_z.setter
    def pos_z(self, value):
        '''
        Sets Z
        '''
        self._pos_z = value

    def position(self):
        '''
        Return the three components of the position
        '''
        return self.pos_x, self.pos_y, self.pos_z

    def dist(self, other):
        '''
        Euclidean distance between two objects
        '''
        pos_x2 = (self.pos_x-other.pos_x)**2.
        pos_y2 = (self.pos_y-other.pos_y)**2.
        pos_z2 = (self.pos_z-other.pos_z)**2.
        return sqrt(pos_x2+pos_y2+pos_z2)

    def transform(self, transformation_matrix):
        '''
        Transforms the position of the point according to the transformation matrix
        '''
        translation = transformation_matrix[:3][3]
        rotation = transformation_matrix[0:3][0:3]
        new_x = translation[0]+rotation[0][0]*self.pos_x+ \
                rotation[0][1]*self.pos_y+rotation[0][2]*self.pos_z
        new_y = translation[1]+rotation[1][0]*self.pos_x+ \
                rotation[1][1]*self.pos_y+rotation[1][2]*self.pos_z
        new_z = translation[2]+rotation[2][0]*self.pos_x+ \
                rotation[2][1]*self.pos_y+rotation[2][2]*self.pos_z
        self.pos_x = new_x
        self.pos_y = new_y
        self.pos_z = new_z
        return self
