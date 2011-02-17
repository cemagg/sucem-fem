import pylab 
from matplotlib.collections import LineCollection
from matplotlib.colors import ColorConverter
colorConverter = ColorConverter()

def draw_2D_edges(edge_coords, axes=None, color='y', show=True):
    if not axes: axes = pylab.axes()
    lc = LineCollection(edge_coords, colors=color, linestyle='solid')
    xc, yc = edge_coords[:,:,0], edge_coords[:,:,1]
    xmin, xmax, ymin, ymax = xc.min(), xc.max(), yc.min(), yc.max()
    xlen = abs(xmax-xmin) ; ylen = abs(ymax-xmin)
    fac = 0.05
    axes.set_xlim((xmin-fac*xlen, xmax+fac*xlen))
    axes.set_ylim((ymin-fac*ylen, ymax+fac*ylen))
    axes.add_collection(lc)
    if show: pylab.show()
    
