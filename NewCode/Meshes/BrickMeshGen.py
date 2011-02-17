import numpy as N

def make_rect_cavity_brick_listmesh(
    a,b,c, edge_len, grid_offset=[0,0,0]):
    cav_dim = N.array((a,b,c), N.float64)
    no_edges = N.ceil(cav_dim/edge_len - 1e-12).astype(N.int32)
    step_size = cav_dim/no_edges
    return {
    'GridDimension' : no_edges+1,
    'GridStepSize' : step_size,
    'GridOffset' : N.array(grid_offset, N.float64),
    }
