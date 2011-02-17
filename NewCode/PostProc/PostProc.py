from itertools import izip
import numpy as N
from numpy import array, zeros, linspace, int32, float64, newaxis

def LocatePoints(mesh, global_coords):
    """
    Return the element number and local coordinate for each point in global_coords
    
    Input Arguments
    ===============
    
    mesh
      Mesh.Mesh or compatible instance
       
    global_coords
      Coordinates to be found, one coordinate triplet per row

    Output 
    ======
    (elnos, local_coords) where

    elnos
      Ordered list of element element numbers of the corresponding global
      coordinate in global_coords

    local_coords
      Ordered list of element local coordinates of the corresponding global
      coordiate in global_coords and element number in elnos
    
    """
    no_coords = len(global_coords)
    elnos = zeros(len(global_coords), int32) - 1
    local_coords = zeros((no_coords, mesh.localCoordLen), float64)
    if hasattr(mesh, "locatePoint"):
        for i,g_coord in enumerate(global_coords):
            elnos[i], local_coords[i] = mesh.locatePoint(global_coords[i])

        return elnos, local_coords

    for i, g_coord in enumerate(global_coords):
        closest_node = mesh.findClosestNode(g_coord)
        node_elnos = mesh.nodeElementConnections[closest_node]
        elno, local_c = LookInElementAndNeigbors(mesh, node_elnos.tolist(), g_coord)
        if elno == -1:      # Sometimes using only the closest node doesn't work
            closest_nodes = mesh.findNodesRadius(g_coord, mesh.maxEdgeLength)
            node_elnos = []
            for node in closest_nodes:
                node_elnos.extend(mesh.nodeElementConnections[node])            
            elno, local_c = LookInElementAndNeigbors(
                mesh, N.unique(node_elnos).tolist(), g_coord)
        elnos[i] = elno
        local_coords[i] = local_c
    return (elnos, local_coords)

def LookInElementAndNeigbors(mesh, elnos, g_coord):
    for els in mesh.elements[elnos].connect2elem:
        elnos.extend(els)
    possible_elements = N.unique(elnos)
    for elno, el in izip(possible_elements, mesh.elements[possible_elements]):
        if el.InElement(g_coord):
            return elno, el.global2local(g_coord)
    return -1, None

def ReconstructPoints(disc_dofs, xyz_coords):
    """
    Reconstruct the field at XYZ coords listed in xyz_coords
    """

    (elnos, local_coords) = LocatePoints(disc_dofs.disc.mesh, xyz_coords)
    return disc_dofs.reconstruct(elnos, local_coords)


def MakeLine(start, stop, no_pts):
    delta = stop - start
    return start + linspace(0, 1, no_pts)[:, newaxis]*delta
