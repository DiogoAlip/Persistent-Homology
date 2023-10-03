import copy
from itertools import combinations

vertices = [[1.0],[0.0],[0,1]]

simplice_0_dim = [0]
simplice_0_dim = [1,2]
simplice_0_dim = [0,1,2]

simplices = [[0],[1],[2],[0,1],[1,2]]

def dimen(simplex):
    return len(simplex)-1

def is_boundary(possible_boundary, simplex):
    if len(possible_boundary) +1 != len(simplex):
        return False
    extra = 0
    for k in range (len(possible_boundary)):
        if possible_boundary[k] == simplex[k+extra]:
            continue
        if extra == 0 and possible_boundary[k] == simplex[k+1]:
            extra = 1
            continue
        return False
    return True

def is_coboundary (possible_coboundary, simplex):
    return is_boundary(simplex, possible_coboundary)

filtracion_por_simplices = [
    [[1, 0], [0, 0], [0, 1]],
    [[0], [1], [2], [0, 1], [1, 2]]
]


#Para el acceso a los elementos de la ﬁltración creamos los siguientes métodos:

def get_vertices(filtration):
    return filtration[0]

def get_vertex(filtration, i):
    return filtration [0] [i]


#Para dibujar complejos simpliciales hasta en 3 dimensiones utilizaremos las siguientes funciones:

def plot_simplicial_complex_0_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_0_dim = simplicial_complex[1][:number_of_vertices]
    return point([get_vertex(simplicial_complex, index[0]) for index in simplex_0_dim], color="red")

def plot_simplicial_complex_1_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_1_dim = []
    for i in range(number_of_vertices, len(simplicial_complex[1])):
        if len(simplicial_complex[1][i]) != 2:
            break
        simplex_1_dim.append(simplicial_complex[1][i])
    return sum(line( [get_vertex(simplicial_complex, index[0]), get_vertex(simplicial_complex, index[1])], color="black"
) for index in simplex_1_dim)

def plot_simplicial_complex_2_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_2_dim = []
    for i in range(number_of_vertices, len(simplicial_complex[1])):
        long = len(simplicial_complex[1][i])
        if long == 3:
            simplex_2_dim.append(simplicial_complex[1][i])
        elif long > 3:
            break
    return sum(polygon([get_vertex(simplicial_complex, index[0]),
                        get_vertex(simplicial_complex, index[1]),
                        get_vertex(simplicial_complex, index[2])])
            for index in simplex_2_dim)

def plot_simplicial_complex_3d (simplicial_complex):
    return sum([plot_simplicial_complex_0_dim(simplicial_complex),
                plot_simplicial_complex_1_dim(simplicial_complex),
                plot_simplicial_complex_2_dim(simplicial_complex)])


#Utilizamos la siguiente función para hacer print de matrices con una ﬁla por línea:

def matrix_print(cadena, M_t):
    print(cadena)
    for i in M_t:
        print(i)


def matrix_print_traspose(cadena, M_t):
    print(cadena)
    for i in traspose(M_t):
        print(i)

def get_boundary_matrix_t(rows, columns):
    return [[1 if is_boundary(rows[j], columns[i]) else 0
                    for j in range(len(rows))]
                    for i in range(len(columns))]

def get_coboundary_matrix_t(rows, columns):
    coboundary_t = []
    for i in range(len(columns)):
        coboundary_t.append([])
        for j in range(len(rows)):
            coboundary_t[i].append(1 if is_coboundary(rows[j], columns[i]) else 0)
    return coboundary_t

def fill_with_zeros(zeros_left, cadena, zeros_right):
    return [0 for _ in range(zeros_left)] + cadena + [0 for _ in range(zeros_right)]

def same_dimension(simplex_1, simplex_2):
    return len(simplex_1) == len(simplex_2)

def boundary_matrix_t(simplicial_complex):
    simplices_list = simplicial_complex[1]
    number_of_vertices = len(simplicial_complex[0])
    number_of_simplices = len(simplices_list)
    
    M_boundary_t = []
    empty_row = [0 for _ in range(number_of_simplices)]
    for _ in range(number_of_vertices):
        M_boundary_t.append(empty_row)

    start_row_index = 0
    start_column_index = number_of_vertices
    for column_index in range(start_column_index, number_of_simplices-1):
        if not same_dimension(simplices_list[column_index], simplices_list[column_index+1]):
            m_boundary_t = get_boundary_matrix_t(simplices_list[start_row_index: start_column_index], simplices_list[start_column_index:column_index+1])

#Fill with 0 boundary matrix
            for col_index in range(len(m_boundary_t)):
                M_boundary_t.append(fill_with_zeros(start_row_index,
                m_boundary_t[col_index], number_of_simplices - start_column_index))
            start_row_index = start_column_index
            start_column_index = column_index+1

    m_boundary_t = get_boundary_matrix_t(simplices_list[start_row_index: start_column_index], simplices_list[start_column_index:])

    #Fill with 0 boundary matrix

    for col_index in range(len(m_boundary_t)):
        M_boundary_t.append(fill_with_zeros(start_row_index, m_boundary_t[col_index], number_of_simplices - start_column_index))
        return M_boundary_t

def low(col):
    for possible_pivot in range(len(col)-1, -1, -1):
        if col[possible_pivot] == 1:
            return possible_pivot
    return -1

def lows(R_t):
    if len(R_t) == 0: return []
    pivots = []
    columns_to_search = [*range(len(R_t))]
    for row_index in range(len(R_t[0])-1, -1, -1):
        for i in range(len(columns_to_search)):
            if R_t[columns_to_search[i]][row_index] != 0:
                pivots.insert(0,row_index)
                columns_to_search.pop(i)
                break
    return pivots

def sum_cols_Z2(first_col, second_col):
    return [1 if first_col[i] != second_col[i] else 0
            for i in range(len(first_col))]

def create_V_matrix(dimen):
    return [[1 if i == j else 0 for j in range(dimen)] for i in range(dimen)]

def get_persistence_pairs_t(R_t):
    persistence_pairs = []
    for j in range(len(R_t)):
        lower = low(R_t[j])
        if lower != -1:
            persistence_pairs.append([lower, j])
    return persistence_pairs

def get_essential_indices_t(R_t, persistence_pairs):
    essential_indices = []
    birth = [pair[0] for pair in persistence_pairs]
    index = 0
    for essential in range(len(R_t)):
        if low(R_t[essential]) == -1 and essential not in birth:
            essential_indices.append(essential)
    return essential_indices

def divide_by_dimension(simplexwise_filtration):
    filtration_by_dimension = []
    p_simplices = []
    for i in range(len(simplexwise_filtration)-1):
        p_simplices.append(simplexwise_filtration[i].copy())
        if dimen(simplexwise_filtration[i]) != dimen(simplexwise_filtration[i+1]):
            filtration_by_dimension.append(p_simplices.copy())
            p_simplices = []
    p_simplices.append(simplexwise_filtration[-1].copy())
    filtration_by_dimension.append(p_simplices.copy())
    return filtration_by_dimension

def persistence_diagram_radious_points(number_of_simplices, diameters, persistence_pairs, essential_indices):
    points = []
    maximum = max(diameters)
    for pair in persistence_pairs:
        points.append([diameters[pair[0]], diameters[pair[1]]])

    for essential_index in essential_indices:
        points.append([diameters[essential_index], maximum])
    
    return points

def dimension_change_indices (split_filtration):
    indices = []
    summand = 0
    for d in split_filtration:
        indices.append(len(d)+summand)
        summand = indices[-1]
    return indices

def simplex_diameter (simplex, distance_m):
    if len(simplex) == 1:
        return 0
    return max([ distance_m[comb[0]][comb[-1]] for comb in combinations(simplex, 2) ])

def plot_persistence_diagram(simplexwise_filtration, persistence_pairs, essential_indices):
    split_filtration = divide_by_dimension(simplexwise_filtration)
    dimension_change_list = dimension_change_indices(split_filtration)
    
    number_of_pp = len(persistence_pairs)
    max_index = len(simplexwise_filtration)-1
    step = max_index/10
    
    points_by_dimension = [[] for _ in range(len(split_filtration))]
    for pair in persistence_pairs:
        for index in range(len(dimension_change_list)):
            if pair[0] < dimension_change_list[index]:
                points_by_dimension[index].append(pair)
                break
    for essential in essential_indices:
        for index in range(len(dimension_change_list)):
            if essential < dimension_change_list[index]:
                points_by_dimension[index].append(
                    [essential, max_index+(step/10)]
                )
                break

    colors = ["red", "blue", "green"]

    for i in range(3, len(points_by_dimension)):
        colors.append((random(),random(),random()))
    return sum(
        point(points_by_dimension[i], color=colors[i], size=15, legend_label=i)
        for i in range(len(points_by_dimension))
    )+line(
        [[0, max_index+(step/10)],[max_index, max_index+(step/10)]],
        color="black", linestyle='--',legend_label='infinity'
    )+line(
        [[0, 0],[max_index, max_index]],
        color="black", linestyle='--',legend_label='x=y'
    )

def count_and_remove_occurrences (element, element_list):
    count = 0
    index = 0
    while index < len(element_list):
        if element_list[index] == element:
            count += 1
            element_list.pop(index)
        else:
            index += 1
    return count

def test ():
    element = 2
    lista = [1,2,3,4,2]
    lista_exp = [1,3,4]
    count_exp = 2
    count = count_and_remove_occurrences (element, lista)
    assert lista == lista_exp
    assert count == count_exp

test()

def test ():
    element = 2
    lista = [1,2,2,4,2]
    lista_exp = [1,4]
    count_exp = 3
    count = count_and_remove_occurrences (element, lista)
    assert lista == lista_exp
    assert count == count_exp
test()

def plot_persistence_diagram_radious(simplexwise_filtration, distance_m, persistence_pairs, essential_indices):
    number_of_simplices = len(simplexwise_filtration)
    
    diameters = [simplex_diameter(simplex, distance_m) for simplex in
    → simplexwise_filtration]
    
    max_diameter = max(diameters)
    step = max_diameter/10

    points = persistence_diagram_radious_points(
        number_of_simplices,
        diameters,
        persistence_pairs,
        essential_indices
    )

    number_of_pp = len(persistence_pairs)
    
    split_filtration = divide_by_dimension(simplexwise_filtration)
    dimension_change_list = dimension_change_indices(split_filtration)

    points_by_dimension = [[] for _ in range(len(split_filtration))]
    point_index = 0
    while point_index < number_of_pp:
        for index in range(len(dimension_change_list)):
            if persistence_pairs[point_index][0] < dimension_change_list[index]:
                points_by_dimension[index].append(points[point_index])
                break
        point_index += 1

    while point_index < len(points):
        for index in range(len(dimension_change_list)):
            if essential_indices[point_index - number_of_pp] <
            → dimension_change_list[index]:
            points_by_dimension[index].
                → append([points[point_index][0],points[point_index][1]+(step/10)])
                break
        point_index += 1

    colors = ["red", "blue", "green"]
    for i in range(3, len(points_by_dimension)):
        colors.append((random(),random(),random()))
    
    points_count = copy.deepcopy(points)
    multiplicity = []
    point_index = 0
    while point_index < len(points_count):
        p = points_count[point_index]
        count = count_and_remove_occurrences (p, points_count)
        if count != 1:
            multiplicity.append([count, [p[0],p[1]+(step/4)]])
    return line(
        [[0, max_diameter+(step/10)],[max_diameter, max_diameter+(step/10)]],
        color="black", linestyle='--',legend_label='infinity'
    )+line(
        [[0, 0],[max_diameter, max_diameter]],
        color="black", linestyle='--',legend_label='x=y'
    )+sum(
        point(points_by_dimension[i], color=colors[i], size=15, legend_label=i)
        for i in range(len(points_by_dimension))
    )+sum(
        text(multi[0], (multi[1]), fontsize='xx-small') for multi in multiplicity
    )

def traspose(matrix):
    if len(matrix) == 0: return []
    return [[matrix[i][j] for i in range(len(matrix))] for j in
    → range(len(matrix[0]))]

def reduce_column_if_necessary(R_t, V_t, i):
    actual_low = low(R_t[i])
    j = 0
    while actual_low != -1 and j < i:
        if low(R_t[j]) == low(R_t[i]):
            R_t[i] = sum_cols_Z2 (R_t[i], R_t[j])
            V_t[i] = sum_cols_Z2 (V_t[i], V_t[j])
            actual_low = low(R_t[i])
            j = 0
        j += 1
    return R_t, V_t

def pHcol_Z2(M_boundary_t):
    assert len(M_boundary_t) > 0
    R_t = copy.deepcopy(M_boundary_t)
    number_of_columns = len(R_t)
    V_t = create_V_matrix(number_of_columns)
    for i in range(number_of_columns):
        R_t, V_t = reduce_column_if_necessary(R_t, V_t, i)
    return R_t, V_t

def indices_with_low_i (R_t, i):
    assert len(R_t) > 0
    indices = []
    for j in range(len(R_t)):
        if low(R_t[j]) == i:
            indices.append(j)
    return indices

def pHrow_Z2(M_boundary_t):
    assert len(M_boundary_t) > 0
    R_t = copy.deepcopy(M_boundary_t)
    dimen = len(R_t)
    V_t = create_V_matrix(dimen)
    for i in range(len(R_t[0])-1, -1, -1):
        indices = indices_with_low_i (R_t, i)
        if len(indices) > 0:
            reduced_col_index = indices[0]
            for j in indices[1:]:
                R_t[j] = sum_cols_Z2 (R_t[j], R_t[reduced_col_index])
                V_t[j] = sum_cols_Z2 (V_t[j], V_t[reduced_col_index])
    return R_t, V_t


cocycle = [0, [[1],[2]]]

cocycles = [[1,[[0,1]]], [1,[[1]]], [0,[[0],[1]]]]

def marked(cocycle):
    return cocycle[0] == 1

def unmarked(cocycle):
    return cocycle[0] == 0

def mark(cocycle):
    cocopy = copy.deepcopy(cocycle)
    cocopy[0] = 1
    return cocopy

def get_chain(cocycle):
    return cocycle[1]

#33 -> 61

def simplex_indices_with_this_coboundary (possible_coboundary, simplices_list):
    boundary_indices = []
    for j in range(len(simplices_list)):
        if is_coboundary(possible_coboundary, simplices_list[j]):
            boundary_indices.append(j)
    return boundary_indices

def is_coboundary_of_the_chain(possible_coboundary, simplices_cochain):
    indices = simplex_indices_with_this_coboundary(possible_coboundary,
    → simplices_cochain)
    return len(indices) % 2 != 0

def cocycles_indices_with_this_coboundary(possible_coboundary, cocycles):
    cocycles_indices = []
    for j in range(len(cocycles)):
        if unmarked(cocycles[j]) and
            → is_coboundary_of_the_chain(possible_coboundary, get_chain(cocycles[j])):
            cocycles_indices.append(j)
    return cocycles_indices

def pCoh_Z2(simplexwise_filtration):
    number_of_simplices = len(simplexwise_filtration)
    
    persistence_pairs=[]
    cocycles = []
    birth = []
    
    for i in range(number_of_simplices):
        new_simplex = simplexwise_filtration[i]
→
        former_cocycles_indices = cocycles_indices_with_this_coboundary(new_simplex, cocycles)

        if len(former_cocycles_indices) == 0:
            cocycles.insert(0,[0,[new_simplex]])
            birth.insert(0,i)
        else:
            p = former_cocycles_indices[0]
            for former_cocycle_index in former_cocycles_indices[1:]:
                cocycles[former_cocycle_index][1] =
                    → cocycles[former_cocycle_index][1] + cocycles[p][1]
            cocycles[p] = mark(cocycles[p])
            persistence_pairs.append([birth[p], i])
            cocycles.insert(0,[1,[new_simplex]])
            birth.insert(0,i)
    return persistence_pairs, cocycles, birth

def get_essential_indices_coh(cocycles, birth):
    essential_indices = []
    for i in range(len(cocycles)):
        if not marked(cocycles[i]):
            essential_indices.append(birth[i])
    return essential_indices

#35 -> 62 page