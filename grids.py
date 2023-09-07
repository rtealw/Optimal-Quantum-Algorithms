import numpy as np
import quantum_query_optimizer as qqo

def get_boolean_labels(nrows, ncols):
    labels = []
    for is_horizontal in [True, False]:
        for row_index in range(nrows):
            for col_index in range(ncols):
                edge_to_start = (row_index == 0 and col_index == 0)
                horizontal_to_first_col = (is_horizontal and col_index == 0)
                vertical_to_first_row = (not is_horizontal and row_index == 0)
                if not (edge_to_start or horizontal_to_first_col or vertical_to_first_row):
                    labels += [(is_horizontal, row_index, col_index)]
    return labels

def is_directed_path(labels, x_str, nrows, ncols):
    x = [int(c) for c in x_str]
    instance = dict(zip(labels, x))
    to_explore = [(0,0)]
    # Depth first search from (0, 0) to (nrows - 1, ncols - 1)
    while len(to_explore) > 0:
        current_vertex = to_explore.pop()
        if current_vertex == (nrows - 1, ncols - 1):
            return '1'
        # Check vertical out
        if current_vertex[0] + 1 < nrows: # Check not top row
            vertical_edge = (False, current_vertex[0]+1, current_vertex[1])
            edge_present = instance[vertical_edge] == 1
            if edge_present:
                next_vertical_vertex = (current_vertex[0]+1, current_vertex[1])
                to_explore += [next_vertical_vertex]

        # Check horizontal out
        if current_vertex[1] + 1 < ncols: # Check not right column
            horizontal_edge = (True, current_vertex[0], current_vertex[1]+1)
            edge_present = instance[horizontal_edge] == 1
            if edge_present:
                next_horizontal_vertex = (current_vertex[0], current_vertex[1]+1)
                to_explore += [next_horizontal_vertex]
    return '0'

nrows = 2
ncols = 3
labels = get_boolean_labels(nrows, ncols)
n = len(labels)
D = qqo.getDAll(n)
E = [is_directed_path(labels, x, nrows, ncols) for x in D]
solutions = qqo.runSDP(D=D, E=E)
