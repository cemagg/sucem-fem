from __future__ import division

import numpy as N
from scipy.io import numpyio
from itertools import izip

def read_filelog(fileobj):
    typecode, rowlen = fileobj.readline().split()
    rowlen = int(rowlen)
    tmparr = []
    while True:
        dline = numpyio.fread(fileobj, rowlen, typecode)
        if len(dline) == 0: break
        if len(dline) < rowlen: raise Exception("Incomplete line read")
        tmparr.append(dline)

    return N.array(tmparr, dtype=N.typeDict[typecode])

def read_filelog_lines(fileobj):
    typecode, rowlen = fileobj.readline().split()
    rowlen = int(rowlen)
    while True:
        dline = numpyio.fread(fileobj, rowlen, typecode)
        if len(dline) == 0: break
        if len(dline) < rowlen: raise Exception("Incomplete line read")
        yield dline


def maxima_mesh(mesh):
    yield ''
    yield 'mesh_verts:rationalize(matrix('
    nodes = iter(mesh.nodes)
    yield str(list(nodes.next()))
    for vert in nodes: yield ',' + str(list(vert))
    yield '))$'
    yield ''
    yield 'tets:['
    tets = iter(mesh.elements)
    yield str(list(tets.next().nodes+1))
    for tet in tets: yield ',' +str(list(tet.nodes+1))
    yield ']$'

def maxima_write_mesh(mesh, f):
    def add_newlines(string_seq):
        for x in string_seq:
            yield x
            yield '\n'
    for x in add_newlines(maxima_mesh(mesh)): f.writelines(x)

def maxima_write_dofs(dofs, f):
    f.writelines(['DOFs: rationalize(', str(dofs.tolist()).replace('j', '*%i'), ')$'])
    f.write('\n')

def maxima_write_dofmap(disc, f):
    """
    Will only yield correct results if there are no constrained DOFs!!!
    """
    f.write('dofmap:[')
    els = iter(disc.elements)
    to_str = lambda el: str((disc.permuter.permuteElement(el)[1]+1).tolist())
    f.write(to_str(els.next()))
    for el in els: f.writelines([',', to_str(el)])
    f.write(']$\n')


def UNV_writemesh(mesh, f):
    # Mesh Nodes
    print >> f, '    -1'
    print >> f, '  2411'
    for i, nc in enumerate(mesh.nodes):
        print >> f, "% 10d% 10d% 10d% 10d" % (i+1,1,1,11)
        print >> f, "  %#- 23.17G  %#- 23.17G  %#- 23.17G" % tuple(list(nc))
    print >> f, '    -1'
    # Mesh Elements
    print >> f, '    -1'
    print >> f, '  2412'
    for i, el in enumerate(mesh.elements): 
        print >> f, "% 10d% 10d% 10d% 10d% 10d% 10d" % (i+1, 111, 1, 1, 7, 4)
        print >> f, "% 10d% 10d% 10d% 10d" % tuple(el.nodes+1)
    print >> f, '    -1'

def mm_write(fd, mat):
    coo_mat = mat.tocoo()
    fd.write('%%MatrixMarket matrix coordinate real general\n')
    fd.write('%d %d %d\n' % (coo_mat.shape[0], coo_mat.shape[1], coo_mat.nnz))
    for row, col, val in izip(coo_mat.row+1, coo_mat.col+1, coo_mat.data):
        fd.write('%d %d %.16e\n' % (row, col, val))
