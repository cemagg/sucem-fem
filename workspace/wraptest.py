import os
from scipy import *
from femfeko import *

def get_dofmap(map):
    n = coupled_td.get_n(map)
    return coupled_td.wrap_dofmap(n, map)
    
femfeko.run_femfeko()

# n = coupled_td.get_n('e')
# print coupled_td.wrap_dofmap(n, 'e')
# n = coupled_td.get_n('d')
# print coupled_td.wrap_dofmap(n, 'd')
# n = coupled_td.get_n('h')
# print coupled_td.wrap_dofmap(n, 'h')
# n = coupled_td.get_n('b')
# print coupled_td.wrap_dofmap(n, 'b')

geomwrap.init_geom()

print geomwrap.vertex_coords
