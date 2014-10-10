'''
Implements functions related to TM align (and eventually TM align itself, hopefully...

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

DIMENSION = 3 # don't like to hard-code

from point import Point
from math import sqrt, atan2, cos, sin

def compute_gl(distances,d0):
    '''
    Computes the GL index given the distances
    '''
    inverse_d0 = 1./d0
    ret = 0
    for dist in distances:
        ret += 1./(1+(dist*inverse_d0)**2)
    return ret

def get_gl(protein1, protein2, alignment, minimum_distance):
    '''
    Gets the GL index between the proteins with the given alignment
    '''
    def get_fit_matrix(coord1, coord2, prev_distances, min_dist):
        indices = []
        while len(indices) < 3 and len(coord1) > 3:
            indices = [i for i in range(len(distances)) \
                    if distances[i] < min_dist]
            md2 = min_dist*min_dist+.5
            min_dist = sqrt(md2)
        new_coord1 = [coord1[x] for x in indices]
        new_coord2 = [coord2[x] for x in indices]
        return fit(new_coord1, new_coord2)

    def get_iteration_gl(orig1, orig2, prev_dists, index_dist, gl_dist):
        if prev_dists == []:
            fit_matrix = fit(orig1, orig2)
        else:
            fit_matrix = get_fit_matrix(orig1, orig2, \
                prev_dists, mini_dist)
        coord1 = [x.transform(fit_matrix) for x in orig1]
        dists = [coord1[i].dist(orig2[i]) \
                for i in range(len(coord1))]
        return compute_gl(dists, minimum_distance), dists


    list_of_aligned = sorted(alignment.keys())
    original_coord1 = [Point((alignment[x].pos_x, alignment[x].pos_y, \
            alignment[x].pos_z)) for x in list_of_aligned]
    original_coord2 = [Point((x.pos_x, x.pos_y, x.pos_z)) \
            for x in list_of_aligned]
    
    #First iteration
    gl1, distances = get_iteration_gl(original_coord1, original_coord2, \
            [], 0, minimum_distance)
    
    #Second iteration
    gl2, distances = get_iteration_gl(original_coord1, original_coord2, \
            distances, minimum_distance, minimum_distance)

    #Third iteration
    md2 = min_dist*min_dist+1.
    gl3 = get_iteration_gl(original_coord1, original_coord2, \
            distances, sqrt(md2), minimum_distance)
    
    return max([gl1, gl2, gl3])

def determinant(matrix):
    '''
    Returns the determinant of the matrix recursively
    '''
    length = len(matrix)
    assert length == len(matrix[0]), "Matrix must be square."
    if length == 1:
        return matrix[0][0]
    else:
        ret = 0
        for j in range(length):
            ret += (-1)**j*matrix[0][j]*determinant(minor(matrix,0,j))
        return ret

def minor(matrix,row,column):
    '''
    Returns the minor to the element i,j; i.e. the matrix without row i and column j
    '''
    length = len(matrix)
    assert length == len(matrix[0]), "Matrix must be square."
    return [[matrix[i][j] for j in range(length) if j!=column] \
            for i in range(length) if i != row]

def matrix_product(matrix, other):
    '''
    Returns the product of the two matrices
    '''
    final_rows = len(matrix)
    final_columns = len(other[0])
    same_dim = len(matrix[0])
    assert same_dim == len(other), "Matrix 1 must have the same number of columns as rows has matrix 2."
    ret = []
    for i in range(final_rows):
        ret.append([])
        for j in range(final_columns):
            ret[i].append(sum([matrix[i][k]*other[k][j] \
                    for k in range(same_dim)]))
    return ret

def matrix_subtract(matrix,other):
    rows = len(matrix)
    cols = len(matrix[0])

    return [[matrix[i][j]-other[i][j] for j in range(cols)] for i in range(rows)]

def tolerate(matrix,tolerance)
    '''
    Repeated procedure in matrix treatment in the original code (not sure why it is like this or what it does exactly)
    '''
    diagonal = sum([matrix[i][1]*matrix[i][0]\
            for i in range(DIMENSION)])
    for i in range(DIMENSION):
        matrix[i][1] -= diagonal_correlation*matrix[i][0]
    local_p = sum([x*x for x in [matrix[y][1] \
            for y in range(DIMENSION)]])
    if local_p <= tolerance:
        helper_m = [matrix[y][0] for y in range(DIMENSION)]]
        local_v = max([abs(x) for x in helper_m])
        local_j = helper_m.index(local_v)
        local_p = sqrt(sum([matrix[k][0]*matrix[k][0] \
                for k in range(DIMENSION) if k != local_j]))
        if local_p <= tolerance: # should stop computing rotation, but I don't know how...
            stop_rotation == True
        matrix[local_j][1] = 0
        matrix[(local_j+1)%DIMENSION][1] = \
                -matrix[(local_j+2)%DIMENSION][0]
        matrix[(local_j+2)%DIMENSION][1] = \
                matrix[(local_j+1)%DIMENSION][0]
    else:
        for i in range(DIMENSION):
            matrix[i][1]/=sqrt(local_p)
    return matrix, stop_rotation



def fit(coord1, coord2):
    '''
    Returns fit_matrix (4x4) that makes the coordinates in coord1 fit on coord2. The matrix contains the rotation matrix (3x3) on fit_matrix[0:3][0:3], a translation vector on fit_matrix[0:3][3] and [0,0,0,1] on fit_matrix[3].
    '''
    tolerance = 1.E-2 #copied from original

    ret = [[0]*(DIMENSION+1) for k in range(DIMENSION+1)]
    for i in range(DIMENSION+1):
        ret[i][i] = 1

#    unique = False #unicity of the superposition
    stop_rotation = False
    length = len(coord1)

    pos1 = []
    pos2 = []
    for i in range(DIMENSION):
        pos1.append([x[i] for x in coord1])
        pos2.append([x[i] for x in coord2])

    sum1 = []
    sum2 = []
    cross = []
    for i in range(DIMENSION):
        sum1.append(sum(pos1[i]))
        sum2.append(sum(pos2[i]))
        cross[i].append([])
        for j in range(3):
            cross[i].append(sum([pos1[i][k]*pos2[j][k] \
                    for k in range(length)]))

    center1 = [x/length for x in sum1]
    center2 = [x/length for x in sum2]
                
    covariance_matrix = [[0]*DIMENSION for k in DIMENSION]
    for i in range(DIMENSION):
        for j in range(DIMENSION):
            covariance_matrix[i][j] = cross[i][j]-sum1[i]*sum2[j]/length

    det_covariance = determinant(covariance_matrix)
    covariance_square = matrix_product(covariance_matrix, covariance_matrix)
    trace_cov_square = sum([covariance_square[i][i] \
            for i in range(DIMENSION)])
    cofactor_cov_square = sum([determinant(minor(covariance_square[i][i])) \
            for i in range(DIMENSION)])
    det_cov_square = det_covariance*det_covariance
    if trace_cov_square > 0: # rotation part is only computed is trace > 0
        # rotation part
        local_h = trace_cov_square*trace_cov_square - cofactor_cov_square
        local_g = (trace_cov_square*cofactor_cov_square-det_cov_square)/2.- \
                trace_cov_square*local_h
        if local_h > 0: #if local_h <= 0, helper_a is identity matrix
            #prepare helper_a
            sqrt_h = sqrt(local_h)
            local_d = atan2(max([local_h*local_h*local_h-local_g*local_g,0]), \
                    -g)/3.
            cos_theta = sqrt_h*cos(local_d)
            sin_theta = sqrt_h*sqrt(3)*sin(local_d)
            helper_e = [trace_cov_square+cos_theta*2, \
                    trace_cov_square-cos_theta+sin_theta, \
                    trace_cov_square-cos_theta-sin_theta]
            if (helper_e[0]-helper_e[1]) > (helper_e[1]-helper_e[2]):
                m1 = 2
                m0 = 0
            else:
                m1 = 0
                m0 = 2
            for l in range(0,DIMENSION,2):
                adjoint_cov_square = [[0]*DIMENSION for k in range(DIMENSION)]
                for i in range(DIMENSION):
                    adjoint_cov_square[i][i] = helper_e[l]
                adjoint_cov_square = matrix_subtract(covariance_square, \
                        adjoint_cov_square)
                diagonal_adjoint_cov_square = [adjoint_cov_square[i][i] \
                        for i in range(DIMENSION)]
                max_diagonal_index = diagonal_adjoint_cov_square.index( \
                        max(diagonal_adjoint_cov_square))
                trace_adj_cov_square = sum([x*x for x in \
                        diagonal_adjoint_cov_square])
                if trace_adj_cov_square > 0:
                    inv_sqrt_trace = 1./sqrt(trace_adj_cov_square)
                else: #if trace_adj_cov_square <= 0, helper_a must not be multiplied
                    inv_sqrt_trace = 1.
                helper_a[l]=[adjoint_cov_square[i][max_diagonal_index] \
                        *inv_sqrt_trace for i in range(DIMENSION)]
                tolerate_a = [helper_a[:][m0], helper_a[:][m1]]
                tolerate_a, stop_rotation = tolerate(tolerate_a, tolerance)
                for i in range[DIMENSION]:
                    helper_a[i][m0] = tolerate_a[i][0]
                    helper_a[i][m1] = tolerate_a[i][1]
            for i in range(DIMENSION):
                helper_a[i][1] = (-1)**(i+1)*det(minor(helper_a,i,1))
        else:
            helper_a = [[0]*DIMENSION for k in range(DIMENSION)]
            for i in range(DIMENSION):
                helper_a[i][i] = 1
        # do helper_b only if we are going to do rotation
        if not stop_rotation:
            helper_b = [[0]*DIMENSION for k in range(DIMENSION)]
            for l in range(DIMENSION-1):
                row_square_sum = 0
                for i in range(DIMENSION):
                    b[i][l] = sum([covariance_matrix[i][k]*helper_a[k][l] \
                            for k in range(DIMENSION)]
                    row_square_sum += b[i][l]*b[i][l]
                    if row_square_sum > 0:
                        inv_sqrt_rqs = 1./sqrt(row_square_sum)
                        for i in range(DIMENSION):
                            b[i][l]*=inv_sqrt_rqs
            tolerate_b = [helper_b[:][0], helper_b[:][1]]
            tolerate_b, stop_rotation = tolerate(tolerate_b, tolerance)
            for i in range[DIMENSION]:
                helper_b[i][0] = tolerate_b[i][0]
                helper_b[i][1] = tolerate_b[i][1]
            for i in range(DIMENSION):
                helper_b[i][2] = (-1)**(i+2)*det(minor(helper_b,i,1))

            # finally, make rotation matrix
            for i in range(DIMENSION):
                for j in range(DIMENSION):
                    ret[i][j] = sum([b[i][k]*a[j][k] for k in range(DIMENSION)])
    # translation part
    for i in range(DIMENSION):
        ret[DIMENSION][i] = center2[i]-sum([ret[i][k]*center1[k] \
                for k in range(DIMENSION)])
    return ret

