from __future__ import division
import numpy as N
import BrickMeshes as BMs

class PyramTestMesh(object):
    no_steps = 3
    test_local_coords =  N.array([(l1, l2, 1-l1-l5, 1-l2-l5, l5)
                                  for l5 in N.linspace(0,99/100, no_steps)
                                  for l2 in N.linspace(0, 1-l5, no_steps)
                                  for l1 in N.linspace(0, 1-l5, no_steps)],
                                 N.float64)

class SixPyram(PyramTestMesh):
    listmesh = {
        'GridDimension' : N.array([2,2,2], N.int32),
        'GridStepSize' : N.array([1,2,3], N.float64),
        'GridOffset' : N.array([0,0,0], N.float64),
        }
    nodes = N.array([
        [0,0,0],[0,0,3],
        [0,2,0],[0,2,3],
        [1,0,0],[1,0,3],
        [1,2,0],[1,2,3],
        [1/2, 1, 3/2]], N.float64)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    pyram_node_offset = 8
    pyram_edge_offset = 12
    nodesIjk = N.array([[i,j,k]         # Only the brick-related nodes are considered
                        for i in range(nnx)
                        for j in range(nny)
                        for k in range(nnz)], N.int32)
    edgeNodes =  N.array([[0,4], [1,5], [2,6], [3,7],
                          [0,2], [1,3], [4,6], [5,7],
                          [0,1], [2,3], [4,5], [6,7], # End of brick edges
                          [0,8], [1,8], [2,8], [3,8],
                          [4,8], [5,8], [6,8], [7,8]], N.int32)
    elementNodes = N.array([[0,1,2,3,8], [4,5,6,7,8],
                            [0,1,4,5,8], [2,3,6,7,8],
                            [0,2,4,6,8], [1,3,5,7,8]], N.int32)
    elementEdges = N.array([[8,4,5,9,12,13,14,15],
                            [10,6,7,11,16,17,18,19],
                            [8,0,1,10,12,13,16,17],
                            [9,2,3,11,14,15,18,19],
                            [4,0,2,6,12,14,16,18],
                            [5,1,3,7,13,15,17,19]], N.int32)
    elementApexfaces = N.array([[8,4,9,5],
                                [10,6,11,7],
                                [8,0,10,1],
                                [9,2,11,3],
                                [4,0,6,2],
                                [5,1,7,3]], N.int32)
    elementBasefaces = N.array([[0], [1], [2], [3], [4], [5]], N.int32)
    elementNodeCoords = N.array([[[ 0  ,  0  ,  0  ], 
                                  [ 0  ,  0  ,  3  ], 
                                  [ 0  ,  2  ,  0  ], 
                                  [ 0  ,  2  ,  3  ], 
                                  [ 1/2,  1  ,  3/2]],
                                 [[ 1  ,  0  ,  0  ], 
                                  [ 1  ,  0  ,  3  ], 
                                  [ 1  ,  2  ,  0  ], 
                                  [ 1  ,  2  ,  3  ], 
                                  [ 1/2,  1  ,  3/2]],
                                 [[ 0  ,  0  ,  0  ], 
                                  [ 0  ,  0  ,  3  ], 
                                  [ 1  ,  0  ,  0  ], 
                                  [ 1  ,  0  ,  3  ], 
                                  [ 1/2,  1  ,  3/2]],
                                 [[ 0  ,  2  ,  0  ], 
                                  [ 0  ,  2  ,  3  ], 
                                  [ 1  ,  2  ,  0  ], 
                                  [ 1  ,  2  ,  3  ], 
                                  [ 1/2,  1  ,  3/2]],
                                 [[ 0  ,  0  ,  0  ], 
                                  [ 0  ,  2  ,  0  ], 
                                  [ 1  ,  0  ,  0  ], 
                                  [ 1  ,  2  ,  0  ], 
                                  [ 1/2,  1  ,  3/2]],
                                 [[ 0  ,  0  ,  3  ], 
                                  [ 0  ,  2  ,  3  ], 
                                  [ 1  ,  0  ,  3  ], 
                                  [ 1  ,  2  ,  3  ], 
                                  [ 1/2,  1  ,  3/2]]], N.float64)
    apexfaceNodes = N.array([[0,4,8], [1,5,8], [2,6,8], [3,7,8],
                             [0,2,8], [1,3,8], [4,6,8], [5,7,8],
                             [0,1,8], [2,3,8], [4,5,8], [6,7,8]], N.int32)
    apexfaceEdges = N.array([[0,16,12], [1,17,13], [2,18,14], [3,19,15],
                             [4,14,12], [5,15,13], [6,18,16], [7,19,17],
                             [8,13,12], [9,15,19], [10,17,16], [11,19,18]],
                            N.int32)
                         
    basefaceNodes = N.array([
        [0,1,3,2], [4,5,7,6],           # x-directed face-normals
        [0,1,5,4], [2,3,7,6],           # y-directed
        [0,2,6,4], [1,3,7,5]            # z-directed
        ], N.int32)
    brickElementEdges = BMs.OneBrick.elementEdges
    brickElementNodes = BMs.OneBrick.elementNodes
    
class TwelvePyram(PyramTestMesh):
    listmesh = {
        'GridDimension' : N.array([2,3,2], N.int32),
        'GridStepSize' : N.array([1,1,1], N.float64),
        'GridOffset' : N.array([0,-1,0], N.float64),
        }
    nodes = N.array([
        [ 0,-1, 0], [ 0,-1, 1],
        [ 0, 0, 0], [ 0, 0, 1],
        [ 0, 1, 0], [ 0, 1, 1],
        [ 1,-1, 0], [ 1,-1, 1],
        [ 1, 0, 0], [ 1, 0, 1],
        [ 1, 1, 0], [ 1, 1, 1],
        [1/2, -1/2, 1/2], [1/2, 1/2, 1/2]], N.float64)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    pyram_node_offset = 12
    pyram_edge_offset = 20
    nodesIjk = N.array([[i,j,k]         # Only the brick-related nodes are considered
                        for i in range(nnx)
                        for j in range(nny)
                        for k in range(nnz)], N.int32)
    edgeNodes =  N.array([
        [0,6], [1,7], [2,8], [3,9], [4,10], [5,11], # x-directed edges
        [0,2], [1,3], [2,4], [3,5], [6,8], [7,9], [8,10], [9,11], # y-directed
        [0,1], [2,3], [4,5], [6,7], [8,9], [10,11],  # z-directed
        [0,12], [1,12], [2,12], [3,12], # Pyramid apex edges
        [6,12], [7,12], [8,12], [9,12],
        [2,13], [3,13], [4,13], [5,13],
        [8,13], [9,13], [10,13],[11,13]], N.int32)
    elementNodes = N.array([[0,1,2,3,12], [6,7,8,9,12],
                            [0,1,6,7,12], [2,3,8,9,12],
                            [0,2,6,8,12], [1,3,7,9,12],
                            [2,3,4,5,13], [8,9,10,11,13],
                            [2,3,8,9,13], [4,5,10,11,13],
                            [2,4,8,10,13], [3,5,9,11,13]], N.int32)
    elementEdges = N.array([[14,6,7,15,20,21,22,23],
                            [17,10,11,18,24,25,26,27],
                            [14,0,1,17,20,21,24,25],
                            [15,2,3,18,22,23,26,27],
                            [6,0,2,10,20,22,24,26],
                            [7,1,3,11,21,23,25,27],
                            [15,8,9,16,28,29,30,31],
                            [18,12,13,19,32,33,34,35],
                            [15,2,3,18,28,29, 32,33],
                            [16,4,5,19,30,31,34,35],
                            [8,2,4,12,28,30,32,34],
                            [9,3,5,13,29,31,33,35]], N.int32)
    elementNodeCoords = N.array([nodes[en] for en in elementNodes], N.float64)
    elementApexfaces = N.array([[8,4,9,5],
                                [10,6,11,7],
                                [8,0,10,1],
                                [9,2,11,3],
                                [4,0,6,2],
                                [5,1,7,3],
                                [20,16,21,17],
                                [22,18,23,19],
                                [20,12,22,13],
                                [21,14,23,15],
                                [16,12,18,14],
                                [17,13,19,15]], N.int32)
    elementBasefaces = N.array([[0],[2],[4],[5],[7],[8],
                                [1],[3],[5],[6],[9],[10]], N.int32)
    apexfaceNodes = N.array([[0,6,12], [1,7,12], [2,8,12], [3,9,12],
                             [0,2,12], [1,3,12], [6,8,12], [7,9,12],
                             [0,1,12], [2,3,12], [6,7,12], [8,9,12],
                             [2,8,13], [3,9,13], [4,10,13], [5,11,13],
                             [2,4,13], [3,5,13], [8,10,13], [9,11,13],
                             [2,3,13], [4,5,13], [8,9,13], [10,11,13]], N.int32)
    apexfaceEdges = N.array([[0,24,20], [1,25,21], [2,26,22], [3,27,23],
                             [6,22,20], [7,23,21], [10,26,24], [11,27,25],
                             [14,21,20], [15,23,22], [17,25,24], [18,27,26],
                             [2,32,28], [3,33,29], [4,34,30], [5,35,31],
                             [8,30,28], [9,31,29], [12,34,32], [13,35,33],
                             [15,29,28], [16,31,35], [18,33,32], [19,35,34]],
                            N.int32)
    basefaceNodes = N.array([
        [0,1,3,2], [2,3,5,4], [6,7,9,8], [8,9,11,10], # x-directed face-normals
        [0,1,7,6], [2,3,9,8], [4,5,11,10],            # y-directed
        [0,2,8,6], [1,3,9,7], [2,4,10,8], [3,5,11,9], # z-directed
        ], N.int32)
    brickElementEdges = BMs.TwoBricks.elementEdges
    brickElementNodes = BMs.TwoBricks.elementNodes
