from __future__ import division
import numpy as N

__all__ = ('OneBrick', 'TwoBricks', 'TwoBricksX', 'FourBricks')


class TestBrick(object):
    no_steps = 3
    test_local_coords =  N.array([(i,j,k,no_steps-i,no_steps-j,no_steps-k)
                                  for i in range(no_steps+1)
                                  for j in range(no_steps+1)
                                  for k in range(no_steps+1)], N.float64)/no_steps

 
class OneBrick(TestBrick):
    listmesh = {
                'GridDimension' : N.array([2,2,2], N.int32),
                'GridStepSize' : N.array([1,2,3], N.float64),
                'GridOffset' : N.array([0,0,0], N.float64),
                }
    nodes = N.array([
        [0,0,0],[0,0,3],
        [0,2,0],[0,2,3],
        [1,0,0],[1,0,3],
        [1,2,0],[1,2,3]], N.float64)
    test_xyz_coords = N.array(
        [TestBrick.test_local_coords[:,0:3]*listmesh['GridStepSize']+listmesh['GridOffset']])
    noGridPoints = len(nodes)
    elementNodes = N.array([[0, 1, 2, 3, 4, 5, 6, 7]], N.int32)
    elementEdges = N.array([[0,1,2,3,4,5,6,7,8,9,10,11]], N.int32)
    elementFaces = N.array([[0,3,1,4,2,5]], N.int32)
    elementConnect2Face = N.array([[-1,-1,-1,-1,-1,-1]], N.int32)
    elementConnect2Elem = N.array([[-1,-1,-1,-1,-1,-1]], N.int32)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    elementGridDimension = N.array([1,1,1], N.int32)
    nodesIjk = N.array([[i,j,k]
                               for i in range(nnx)
                               for j in range(nny)
                               for k in range(nnz)], N.int32)
    edgeNodes =  N.array([[0,4], [1,5], [2,6], [3,7],
                          [0,2], [1,3], [4,6], [5,7],
                          [0,1], [2,3], [4,5], [6,7]], N.int32)
    edgeConnect2Elem = N.array([
        [-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],
        [-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],
        [-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1],[-1,-1,-1,-1]], N.int32)
    edgeGridDim = (4,4,4)               # no of x,y,z directed edges
    faceNodes = N.array([
        [0,1,3,2], [4,5,7,6],           # x-directed face-normals
        [0,1,5,4], [2,3,7,6],           # y-directed
        [0,2,6,4], [1,3,7,5]            # z-directed
        ], N.int32)
    faceEdges = N.array([
        [8,5,9,4], [10,7,11,6],
        [8,1,10,0],[9,3,11,2],
        [4,2,6,0], [5,3,7,1]], N.int32)
    faceConnect2Elem = N.array(
        [[0,-1], [0,-1],
         [0,-1], [0,-1],
         [0,-1], [0,-1]], N.int32)
    faceGridDim = (2,2,2)
    boundaryEdgeSet = set(range(len(edgeNodes)))
    boundaryFaceSet = set(range(len(faceNodes)))
    boundaryFaces = N.array([True]*6, N.bool)
    boundaryEdges = N.array([True]*12, N.bool)

        
class TwoBricks(TestBrick):
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
        [ 1, 1, 0], [ 1, 1, 1]], N.float64)
    test_xyz_coords = N.array(
        [(TestBrick.test_local_coords[:,0:3]*listmesh['GridStepSize']
          + P0) for P0 in nodes[[0,2]]])
    noGridPoints = len(nodes)
    elementNodes = N.array([[0, 1, 2, 3, 6, 7, 8, 9],
                            [2, 3, 4, 5, 8, 9, 10, 11]], N.int32)
    elementEdges = N.array([[0,1,2,3,6,7,10,11,14,15,17,18],
                            [2,3,4,5,8,9,12,13,15,16,18,19]], N.int32)
    elementFaces = N.array([[0,4,7,2,5,8], [1,5,9,3,6,10]], N.int32) 
    elementConnect2Face = N.array([[-1,-1,-1,-1,1,-1],
                                   [-1,4,-1,-1,-1,-1]], N.int32)
    elementConnect2Elem = N.array([[-1,-1,-1,-1,1,-1],
                                   [-1,0,-1,-1,-1,-1]], N.int32)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    elementGridDimension = N.array([1,2,1], N.int32)
    nodesIjk = N.array([[i,j,k]
                               for i in range(nnx)
                               for j in range(nny)
                               for k in range(nnz)], N.int32)
    edgeNodes =  N.array([
        [0,6], [1,7], [2,8], [3,9], [4,10], [5,11], # x-directed edges
        [0,2], [1,3], [2,4], [3,5], [6,8], [7,9], [8,10], [9,11], # y-directed
        [0,1], [2,3], [4,5], [6,7], [8,9], [10,11]  # z-directed
        ], N.int32)
    edgeConnect2Elem = N.array([
        [0,-1,-1,-1],[0,-1,-1,-1], [0,1,-1,-1],
        [0,1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,-1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,-1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,1,-1,-1], [1,-1,-1,-1]], N.int32)
    edgeGridDim = (6,8,6)               # no of x,y,z directed edges
    faceNodes = N.array([
        [0,1,3,2], [2,3,5,4], [6,7,9,8], [8,9,11,10], # x-directed face-normals
        [0,1,7,6], [2,3,9,8], [4,5,11,10],            # y-directed
        [0,2,8,6], [1,3,9,7], [2,4,10,8], [3,5,11,9], # z-directed
        ], N.int32)
    faceEdges = N.array([
        [14,7,15,6], [15,9,16,8], [17,11,18,10], [18,13,19,12],
        [14,1,17,0], [15,3,18,2], [16,5,19,4],
        [6,2,10,0], [7,3,11,1], [8,4,12,2], [9,5,13,3]], N.int32)
    faceConnect2Elem = N.array(
        [[0,-1], [1,-1], [0,-1], [1,-1],
         [0,-1], [0,1], [1,-1],
         [0,-1], [0,-1], [1,-1], [1,-1]], N.int32)
    faceGridDim = (4,3,4)
    boundaryEdgeSet = set(range(len(edgeNodes)))
    boundaryFaceSet = set(range(len(faceNodes))) - set([5])
    boundaryFaces = N.array([i in boundaryFaceSet for i in range(len(faceNodes))],
                            N.bool)
    boundaryEdges = N.array([i in boundaryEdgeSet for i in range(len(edgeNodes))],
                            N.bool)

class TwoBricksX(TestBrick):
    listmesh = {
                'GridDimension' : N.array([3,2,2], N.int32),
                'GridStepSize' : N.array([1,1,1], N.float64),
                'GridOffset' : N.array([0,0,0], N.float64),
                }
    nodes = N.array([
        [0,0,0], [0,0,1],
        [0,1,0], [0,1,1],
        [1,0,0], [1,0,1],
        [1,1,0], [1,1,1],
        [2,0,0], [2,0,1],
        [2,1,0], [2,1,1]], N.float64)
    test_xyz_coords = N.array(
        [(TestBrick.test_local_coords[:,0:3]*listmesh['GridStepSize']
          + P0) for P0 in nodes[[0,4]]])
    noGridPoints = len(nodes)
    elementNodes = N.array([[0, 1, 2, 3, 4, 5, 6, 7],
                            [4, 5, 6, 7, 8, 9,10,11]], N.int32)
    elementEdges = N.array([[0,1,2,3,8,9,10,11,14,15,16,17],
                            [4,5,6,7,10,11,12,13,16,17,18,19]], N.int32)
    elementFaces = N.array([[0,3,7,1,4,8], [1,5,9,2,6,10]], N.int32) 
    elementConnect2Face = N.array([[-1,-1,-1,0,-1,-1],
                                   [3,-1,-1,-1,-1,-1]], N.int32)
    elementConnect2Elem = N.array([[-1,-1,-1,1,-1,-1],
                                   [0,-1,-1,-1,-1,-1]], N.int32)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    elementGridDimension = N.array([2,1,1], N.int32)
    nodesIjk = N.array([[i,j,k]
                               for i in range(nnx)
                               for j in range(nny)
                               for k in range(nnz)], N.int32)
    edgeNodes =  N.array([
        [0,4], [1,5], [2,6], [3,7], [4,8], [5,9], [6,10],[7,11], # x-directed edges
        [0,2], [1,3], [4,6], [5,7], [8,10], [9,11],  # y-directed
        [0,1], [2,3], [4,5], [6,7], [8,9], [10,11] # z-directed
        ], N.int32)
    edgeConnect2Elem = N.array([
        [0,-1,-1,-1], [0,-1,-1,-1], [0,-1,-1,-1], [0,-1,-1,-1],
        [1,-1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,-1,-1,-1], [0,1,-1,-1],
        [0,1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1],
        [0,-1,-1,-1], [0,-1,-1,-1], [0,1,-1,-1],
        [0,1,-1,-1], [1,-1,-1,-1], [1,-1,-1,-1]], N.int32) 
    edgeGridDim = (8,6,6)               # no of x,y,z directed edges
    faceNodes = N.array([
        [0,1,3,2], [4,5,7,6], [8,9,11,10],            # x-directed face-normals
        [0,1,5,4], [2,3,7,6], [4,5,9,8], [6,7,11,10], # y-directed
        [0,2,6,4], [1,3,7,5], [4,6,10,8], [5,7,11,9], # z-directed
        ], N.int32)
    faceEdges = N.array([
        [14,9,15,8], [16,11,17,10], [18,13,19,12], 
        [14,1,16,0], [15,3,17,2], [16,5,18,4], [17,7,19,6],
        [8,2,10,0], [9,3,11,1], [10,6,12,4], [11,7,13,5]], N.int32) 
    faceConnect2Elem = N.array([
        [0,-1], [0,1], [1,-1],
        [0,-1], [0,-1], [1,-1], [1,-1],
        [0,-1], [0,-1], [1,-1], [1,-1]], N.int32)
    faceGridDim = (3,4,4)
    boundaryEdgeSet = set(range(len(edgeNodes)))
    boundaryFaceSet = set(range(len(faceNodes))) - set([1])
    boundaryFaces = N.array([i in boundaryFaceSet for i in range(len(faceNodes))],
                            N.bool)
    boundaryEdges = N.array([i in boundaryEdgeSet for i in range(len(edgeNodes))],
                            N.bool)

class FourBricks(TestBrick):
    listmesh = {
                'GridDimension' : N.array([2,3,3], N.int32),
                'GridStepSize' : N.array([1,2,3], N.float64),
                'GridOffset' : N.array([0,0,0], N.float64),
                }
    nodes = N.array([
        [0,0,0], [0,0,3], [0,0,6],
        [0,2,0], [0,2,3], [0,2,6],
        [0,4,0], [0,4,3], [0,4,6],
        [1,0,0], [1,0,3], [1,0,6],
        [1,2,0], [1,2,3], [1,2,6],
        [1,4,0], [1,4,3], [1,4,6]], N.float64)
    test_xyz_coords = N.array(
        [(TestBrick.test_local_coords[:,0:3]*listmesh['GridStepSize']
          + P0) for P0 in nodes[[0,1,3,4]]])
    noGridPoints = len(nodes)
    elementNodes = N.array([[0,1,3,4,9,10,12,13],
                            [1,2,4,5,10,11,13,14],
                            [3,4,6,7,12,13,15,16],
                            [4,5,7,8,13,14,16,17]], N.int32)
    elementEdges = N.array([[0,1,3,4,9,10,15,16,21,23,27,29],
                            [1,2,4,5,10,11,16,17,22,24,28,30],
                            [3,4,6,7,12,13,18,19,23,25,29,31],
                            [4,5,7,8,13,14,19,20,24,26,30,32]], N.int32)
    elementFaces = N.array([[0,8,14,4,10,15],
                            [1,9,15,5,11,16],
                            [2,10,17,6,12,18],
                            [3,11,18,7,13,19]], N.int32) 
    elementConnect2Face = N.array([[-1,-1,-1,-1,1,2],
                                   [-1,-1,5,-1,1,-1],
                                   [-1,4,-1,-1,-1,2],
                                   [-1,4,5,-1,-1,-1]],N.int32)
    elementConnect2Elem = N.array([[-1,-1,-1,-1,2,1],
                                   [-1,-1,0,-1,3,-1],
                                   [-1,0,-1,-1,-1,3],
                                   [-1,1,2,-1,-1,-1]], N.int32)
    nnx, nny, nnz = listmesh['GridDimension']
    gridDimension = listmesh['GridDimension']
    elementGridDimension = N.array([1,2,2], N.int32)
    nodesIjk = N.array([[i,j,k]
                               for i in range(nnx)
                               for j in range(nny)
                               for k in range(nnz)], N.int32)
    edgeNodes = N.array([
        [0,9], [1,10], [2,11], [3,12], [4,13], # x-directed
        [5,14], [6,15], [7,16], [8,17],        # x
        [0,3], [1,4], [2,5], [3,6], [4,7], [5,8],            # y-directed
        [9,12], [10,13], [11,14], [12,15], [13,16], [14,17], # y
        [0,1], [1,2], [3,4], [4,5], [6,7], [7,8], [9,10],    # z-directed
        [10, 11], [12,13], [13,14], [15,16], [16,17]],       # z
                        N.int32)
    edgeConnect2Elem = N.array([
        [0,-1,-1,-1],[0,1,-1,-1],[1,-1,-1,-1],[0,2,-1,-1],[0,1,2,3], #x
        [1,3,-1,-1], [2,-1,-1,-1],[2,3,-1,-1],[3,-1,-1,-1],
        [0,-1,-1,-1],[0,1,-1,-1],[1,-1,-1,-1], # y
        [2,-1,-1,-1],[2,3,-1,-1],[3,-1,-1,-1],
        [0,-1,-1,-1],[0,1,-1,-1],[1,-1,-1,-1], 
        [2,-1,-1,-1],[2,3,-1,-1],[3,-1,-1,-1],
        [0,-1,-1,-1],[1,-1,-1,-1],      # z
        [0,2,-1,-1],[1,3,-1,-1],
        [2,-1,-1,-1],[3,-1,-1,-1],
        [0,-1,-1,-1],[1,-1,-1,-1],      
        [0,2,-1,-1],[1,3,-1,-1],
        [2,-1,-1,-1],[3,-1,-1,-1]], N.int32) 
    edgeGridDim = (9,12,12)             # no of x,y,z directed edges
    faceNodes = N.array([
        [0,1,4,3], [1,2,5,4], [3,4,7,6], [4,5,8,7],                # x-normal
        [9,10,13,12], [10,11,14,13], [12,13,16,15], [13,14,17,16], # x
        [0,1,10,9], [1,2,11,10], [3,4,13,12],  # y-normal
        [4,5,14,13], [6,7,16,15], [7,8,17,16], # y
        [0,3,12,9], [1,4,13,10], [2,5,14,11],   # z-normal
        [3,6,15,12], [4,7,16,13], [5,8,17,14]], # z
                        N.int32)
    faceEdges = N.array([
        [21,10,23,9], [22,11,24,10], [23,13,25,12], [24,14,26,13],
        [27,16,29,15], [28,17,30,16], [29,19,31,18], [30,20,32,19],
        [21,1,27,0], [22,2,28,1], [23,4,29,3],
        [24,5,30,4], [25,7,31,6], [26,8,32,7],
        [9,3,15,0],[10,4,16,1],[11,5,17,2],
        [12,6,18,3], [13,7,19,4], [14,8,20,5]], N.int32)
    faceConnect2Elem = N.array([
        [0,-1], [1,-1], [2,-1], [3,-1],
        [0,-1], [1,-1], [2,-1], [3,-1],
        [0,-1], [1,-1], [0,2], [1,3], [2,-1], [3,-1],
        [0,-1], [0,1], [1,-1], [2,-1], [2,3], [3,-1]], N.int32)
    faceGridDim = (8, 6, 6)
    boundaryFaceSet = set(range(20)) - set([10,11,15,18])
    boundaryFaces = N.array([i in boundaryFaceSet for i in range(len(faceNodes))],
                            N.bool)
    boundaryEdgeSet = set(range(33)) - set([4])
    boundaryEdges = N.array([i in boundaryEdgeSet for i in range(len(edgeNodes))],
                            N.bool)


