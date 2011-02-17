def constrained_on_boundary (entity):
    """
    Freedom function for PEC boundaries. False for entities on mesh boundary
    """
    # Equivalent to not entity.onBoundary for scalars, and
    #
    # [not onbdry for onbdry in entity.onBoundary]
    #
    # if we're dealing with an array. This allows the same code to work for both.
    return entity.onBoundary < True

E_PECboundary_free_edge = constrained_on_boundary
# Well, the code is the same :)
B_PECboundary_free_face = constrained_on_boundary
allfree = lambda x: True
allconstrained = lambda x: False
