## Copyright (C) 2011 Stellenbosch University
##
## This file is part of SUCEM.
##
## SUCEM is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## SUCEM is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with SUCEM. If not, see <http:##www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
#!/usr/bin/python

# mesh2femmesh.py
# simple script to convert I-DEAS universal (UNV) format
# and the Gmsh .msh format to the internal
# femfeko .femmesh mesh file format
# written and maintained by CEMAGG @US
#
# Authors
# Neilen Marais <nmarais@gmail.com>
# Julian P. Swartz <swartzjp@gmail.com>

import sys
from Scientific.IO.FortranFormat import *

class UNVerror(Exception):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return repr(self.message)


class Element(object):
    shapeid = 'unknown'                 # Geometric shape ID, eg. triangle
    def __init__(self, nodes, material=None):
        self.material=material
        self.nodes = nodes

class Triangle(Element):
    shapeid = 'triangle'

class Tetrahedron(Element):
    shapeid = 'tetrahedron'

class Quadrilateral(Element):
    shapeid = 'quadrilateral'

class Brick(Element):
    shapeid = 'brick'

# 2nd order tetrahedron
class Tetrahedron2(Element):
    shapeid = 'tetrahedron2'
    pass


# JPS: allow the input and output meshfilenames to be passed
# as command line parameters.
# program is then used as:
# $python mesh2femmesh infile outfile

inmeshfilename=sys.argv[1]
outmeshfilename=sys.argv[2]

class InputUNVMesh(object):
    beamElementTypeNos = set([11, 21, 22, 23, 24])
    skipElementTypeNos = set([41])
    def __init__(self, inmeshfile, outmesh):
        self.inmeshfile = inmeshfile
        self.outmesh = outmesh
        # Map of UNV file node labels to FEMFEKO global node numbers
        self.__node_map = {}
        # Map of UNV global element numbers to FEMFEKO face and element numbers
        self.__element_number_map = {
            Triangle: {},
            Tetrahedron: {},
            Quadrilateral: {},
            Brick: {},
            }
        self.__UNV_blocks = {
            '2411' : self.__parse_Nodes_Block,
            '2412' : self.__parse_Elements_Block,
            '2477' : self.__parse_Element_Properties_Block
            }
        self.__UNV_element_map = {
            44: Quadrilateral,
            91: Triangle,
            111: Tetrahedron,
            115: Brick,
            }
        pass

    def __block_Generator(self):
        """
        Returns a generator for blocks. Will loop over a UNV file block,
        yielding lines, stopping at the end of a block. inmeshfile needs
        to be positioned on the line following the block start marker.
        """
        for line in self.inmeshfile:
            if line.strip() == '-1' :       # End of block marker
                return                      
            yield line
            pass
        pass

    def __parse_Element_Properties_Block(self):
        """
        Parse element material properties block
        """
        print 'block_materials'
        
        for line in self.__block:
            # Header containing the material index an no of elements affected
            block_header = line.split()
            material_index = int(block_header[0])
            # list of elements that belong to material_index
            element_list = []
            # Number of elements listed for material_index
            no_elements = int(block_header[7])
            # Odd number of elements will require special treatment of last line
            if no_elements % 2 == 0: no_elements_odd = False
            else: no_elements_odd = True
            # Two elements listed per line. Note integer division returns
            # the floor of the result, so the last line will be unread in the
            # case of an odd number of elements
            no_lines = no_elements/2    
            self.__block.next()         # Skip material group name

            for x in xrange(no_lines):
                split_line = self.__block.next().split()
                elementA = int(split_line[1]) 
                elementB = int(split_line[5]) 
                element_list.extend([elementA, elementB])
                pass

            # Read the last element for the case of an odd number
            if no_elements_odd:
                elementA = self.__block.next().split()[1] 
                element_list.extend([elementA])
                pass

            # The element list contains UNV global element numbers,
            # and no information of the element type. The next for
            # loop, loops over the known element_type s (as found in
            # self.__element_number_map), and finds the FEMFEKO
            # equivalent numbers.
            for element_type, no_map in self.__element_number_map.iteritems():
                # no_map = UNV <-> FEMFEKO number map for element_type
                shapeid = element_type.shapeid

                # Looks up FEMFEKO element number for elements in
                # element_list. If it is element_type, add it to the
                # corresponding mapping list
                material_elements = [
                    no_map[i]
                    for i in element_list 
                    if no_map.has_key(i)
                    ]
                # Set the material index of the elements in material_elements
                self.outmesh.setElementProperties(shapeid, material_elements,
                                                  material_index)
                pass
            pass    
        pass


    def __parse_Elements_Block(self):
        """
        Parse elements block for supported element types. The elements
        are added to self.outmesh using the addElement method. All material
        types are initialised to 0.

        Elements are renumbered so that there each element type is
        numbered from 1.
        """
        
        print 'block_elements'

        element_count = {}              # Element counter
        # Initialise element counter for each type to 0
        for element_class in self.__element_number_map:
            element_count[element_class] = 0

        for line in self.__block:
            # 1st element line
            (element_label, element_type) = \
                            [int(i) for i in line.split()][0:2]
            line = self.__block.next()
            # Beam elements have three lines instead of two but are unsupported
            if element_type in self.beamElementTypeNos:
                self.__block.next()
                continue
            element_nodes = [int(i) for i in line.split()] # 2nd element line
            # Renumber nodes from UNV -> FEMFEKO numbers
            element_nodes = [self.__node_map[i] for i in element_nodes]

            try:
                # Create element with material type defaulted to 0
                element = self.__UNV_element_map[element_type](element_nodes, 0)
                element_class = element.__class__
            except KeyError:
                if element_type in self.skipElementTypeNos: continue
                else:
                    raise UNVerror("Unknown element type %i while reading element block"
                                   % element_type)

            # Add element to mesh
            self.outmesh.addElement(element) 
            # Increment the count for this element type
            element_count[element_class] += 1
            # Get UNV <-> FEMFEKO element number map for this element type
            element_map = self.__element_number_map[element_class]
            # Assign mapping between UNV element_label and FEMFEKO element number
            element_map[element_label] = element_count[element_class] - 1
            pass
        
        pass
    
    def __parse_Nodes_Block(self):
        print 'block_nodes'
        # Contains node lable no, and other unused info
        lineA_format = FortranFormat('4I10')
        # The x, y an z coordinates of the node
        lineB_format = FortranFormat('3D25.16')
        node_count = 0
        
        for line in self.__block:
            dataA = FortranLine(line, lineA_format)
            line = self.__block.next()  # Get the next line
            dataB = FortranLine(line, lineB_format)
            self.outmesh.addNode(dataB[:]) # Add node to the mesh

            node_label = dataA[0]       # First int on lineA is node label
            node_count += 1
            self.__node_map[node_label] = node_count-1 # Mapping for this node
            pass
        pass

    def __block_handler(self):
        # The next line is the block type
        blocktype = self.inmeshfile.next().strip()
        # Initialise the block iterator, assuming self.inmesfile is positioned
        # inside a block
        self.__block = self.__block_Generator()    
        if blocktype in self.__UNV_blocks:    # If blocktype known, call function
                                              # from hash 
            self.__UNV_blocks[blocktype]()
        else:                               
            self.__endBlock()                 # Skip to end of block
            pass
    
    def __endBlock(self):
        for line in self.__block:
            pass
        pass


    def readMesh(self):
        # Loop over lines, looking for known UNV data sets
        for line in self.inmeshfile:
            if line.strip() == '-1':            # Beginning of a block
                self.__block_handler()            # Determine block type, and handle it
                pass
            pass
        pass

    def copyNodeMap(self):
        return self.__node_map.copy()
    
    pass

#####
# Start JPS edit
# class to read in Gmsh .msh files
#####

class InputMesh(object):
    def __init__(self, inmeshfile, outmesh):
        self.setBlockTypes()
        self.setElementTypeHandler()
        self.inmeshfile = inmeshfile
        self.outmesh = outmesh
        # Map of Gmsh file node labels to FEMFEKO global node numbers
        self.node_map = {}

    def readMesh(self):
        # Loop over lines, looking for known data sets
        for line in self.inmeshfile:
            if self.startOfBlock(line):
                # Beginning of a block
                print line
                self.block_handler(line.strip()) # Determine block type, and handle it

    def copyNodeMap(self):
        return self.node_map.copy()

    def endBlock(self):
        for line in self.block:
            pass

    def block_Generator(self):
        """
        Returns a generator for blocks. Will loop over a mesh file block,
        yielding lines, stopping at the end of a block. inmeshfile needs
        to be positioned on the line following the block start marker.
        """
        print "in block generator"
        for line in self.inmeshfile:
            if self.endOfBlock(line):       # End of block marker
                return         
            yield line

    def block_handler(self, blocktype):
        # block type is read in previously as the gmsh block delimiter
        print "block type is ", blocktype
        # Initialise the block iterator, assuming self.inmesfile is positioned
        # inside a block
        self.block = self.block_Generator()
        if blocktype in self.blockTypes:    # If blocktype known, call function
                                              # from hash 
            self.blockTypes[blocktype]()
        else:                               
            self.endBlock()                 # Skip to end of block

class InputGmshMesh(InputMesh):
    '''Class to read meshes created by Gmsh. The node and element data is stored in dictionaries
    and written to an output file in the FEMFEKO .femmesh format.'''

    # note that the Gmsh input routines currently only support first and second order tetrahedral elements.
    # the code can however easily be extended to provide this functionality and additional features will be added
    # when required.
    
    # Map of Gmsh global element numbers to FEMFEKO face and element numbers
    element_number_map = {
        Triangle: {},
        Tetrahedron: {},
        Tetrahedron2: {}
        }
    elementTypeMap = {
        2:  Triangle,
        4:  Tetrahedron,
        11: Tetrahedron2  # second order tet
        }

    blockEnd = '$END'
    blockStart = '$'                    # Prolly bit of a hack, but seems to work

    def startOfBlock(self,line):
        return ((line[0] == self.blockStart) and (line[0:4] != self.blockEnd))

    def endOfBlock(self, line):
        return (line[0:4] == self.blockEnd)

    def setBlockTypes(self):
        self.blockTypes = {
            '$NOD' : self.parse_Nodes_Block,
            '$ELM' : self.parse_Elements_Block
            }

    def setElementTypeHandler(self):
        self.elementTypeHandler = {
            4: self.gmshTet,
            11: self.gmshTet2
            }

    def gmshTet(self, splitLine):
        # read node list for first order element
        element_material = int(splitLine[0])
        element_nodes = [int(i) for i in splitLine[3:7]] # element node numbers
        return element_material, element_nodes

    def gmshTet2(self, splitLine):
        # read node list for first order element
        element_material = int(splitLine[0])
        # read node list for second order element
        element_nodes = [int(i) for i in splitLine][3:13] # element node numbers
        # the midpoint nodes must be re-ordered to match the order convention of FEMFEKO
        temp_nodes = element_nodes[:]
        
        element_nodes[5] = temp_nodes[6]
        element_nodes[6] = temp_nodes[7]
        element_nodes[7] = temp_nodes[5]
        element_nodes[8] = temp_nodes[9]
        element_nodes[9] = temp_nodes[8]
        return element_material, element_nodes
    
    def parse_Elements_Block(self):
        """
        Parse elements block for supported element types. The elements
        are added to self.outmesh using the addElement method. All material
        types are initialised to 0.

        Elements are renumbered so that there each element type is
        numbered from 1.
        """
        
        print 'block_elements'

        line = self.block.next()  # skip number of elements line

        element_count = {}              # Element counter
        # Initialise element counter for each type to 0
        for element_class in self.element_number_map:
            element_count[element_class] = 0

        for line in self.block:
            #print line.strip() # the strip() removes the newline character
            #  element line
            #temp = line.split()
            (element_label, element_type) = \
                            [int(i) for i in line.split()][0:2]

            try: element_material, element_nodes = \
                 self.elementTypeHandler[element_type](line.split()[2:])
            except KeyError:
                print "Skipping element type ", element_type
                continue
                
            #print element_nodes
            # Renumber nodes from Gmsh -> FEMFEKO numbers
            element_nodes = [self.node_map[i] for i in element_nodes]
            
            try:
                # Create element with material type defaulted to 0
                element = self.elementTypeMap[element_type](element_nodes, element_material)
                element_class = element.__class__
            except KeyError:
                raise UNVerror("Unknown element type %i while reading element block"
                               % element_type)
            
            # Add element to mesh
            self.outmesh.addElement(element)
            # Increment the count for this element type
            element_count[element_class] += 1
            # Get Gmsh <-> FEMFEKO element number map for this element type
            element_map = self.element_number_map[element_class]
            # Assign mapping between UNV element_label and FEMFEKO element number
            element_map[element_label] = element_count[element_class] - 1
    
    def parse_Nodes_Block(self):
        print 'block_nodes'
        node_count = 0
        line = self.block.next()  # skip number of nodes line
        
        for line in self.block:
            #print line   # feedback for debugging
            temp = line.split()
            dataA = int(temp[0])
            dataB = [float(temp[1]), float(temp[2]), float(temp[3])]
            self.outmesh.addNode(dataB[:]) # Add node to the mesh

            node_label = dataA       #  node label
            node_count += 1
            self.node_map[node_label] = node_count-1 # Mapping for this node

#####
# End JPS edit
#####

class InputGmsh2Mesh(InputGmshMesh):
    blockEnd = '$End'

    def setBlockTypes(self):
        self.blockTypes = {
            '$Nodes' : self.parse_Nodes_Block,
            '$Elements' : self.parse_Elements_Block
            }

    def gmshTet(self, splitLine):
        # read node list for first order element
        no_tags = int(splitLine[0])
        element_material = splitLine[1]
        elNodes = splitLine[no_tags+1:]
        element_nodes = [int(i) for i in elNodes] # element node numbers
        return element_material, element_nodes

    def gmshTet2(self, splitLine):
        raise NotImplemented(
            "Compare to JPS's gmsh v1.0 format code, I (Neilen) am not inclined to test this")

    def setElementTypeHandler(self):
        self.elementTypeHandler = {
            2: self.gmshTet,
            4: self.gmshTet,
            11: self.gmshTet2
            }


class InputNetgenMesh(InputMesh):
    def setBlockTypes(self):
        self.blockTypes = {
            'volumeelements': self.parseElementsBlock,
            'points' : self.parseNodesBlock
            }

    def setElementTypeHandler(self): pass

    def readMesh(self):
        # Loop over lines, looking for known data sets
        for line in self.inmeshfile:
            if self.startOfBlock(line):
                # Beginning of a block
                print line
                # Slight hack, since node coords are specified last in the file
                if line.strip() == 'points':
                    self.block_handler(line.strip()) # Determine block type, and handle it
                    break

        self.inmeshfile.seek(0)
        for line in self.inmeshfile:
            if self.startOfBlock(line):
                # Beginning of a block
                print line
                # Slight hack, since node coords are specified last in the file
                if line.strip() == 'volumeelements':
                    self.block_handler(line.strip()) # Determine block type, and handle it
                    break

        


    def parseNodesBlock(self):
        print 'block Nodes'
        node_count = 0
        line = self.block.next()        # skip number of nodes line
        for line in self.block:
            coords = [float(x) for x in line.split()]
            self.outmesh.addNode(coords)
            node_count += 1
            node_label = node_count
            self.node_map[node_label] = node_count - 1

    def parseElementsBlock(self):
        line = self.block.next()        # skip numer of elements
        element_count = 0
        for line in self.block:
            splitLine = line.split()
            element_count += 1
            element_label = element_count
            element_material = int(splitLine[0])
            element_nodes = [int(f) for f in splitLine[2:6]]
            element_nodes = [self.node_map[i] for i in element_nodes]
            element = Tetrahedron(element_nodes, element_material)
            self.outmesh.addElement(element)
            
    def startOfBlock(self,line):
        return (line.strip() in self.blockTypes)

    def endOfBlock(self, line):
        return (line.strip() == '')

class OutputMesh(object):
    def __init__(self, outmeshfile, inmesh):
        self.outmeshfile = outmeshfile
        self.inmesh = inmesh
        pass

    def writeMesh(self):
        """
        Write inmesh to a FEMFEKO .femmesh format file
        """

        self.writeNodes()
        # write out triangles if any
        if len(self.inmesh.triangles)>0:
            self.writeTris()
        # write out first order tets if any
        if len(self.inmesh.tetraheders)>0:
            self.writeTets()
        # write out second order tets if any
        if len(self.inmesh.tetraheders2)>0:
            self.writeTet2s()
        pass

    def writeTris(self):
        """
        Writes triangles to the femmesh output file
        """

        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK tris"
        print >> outmeshfile, len(inmesh.triangles)

        # first write out tris
        tri_index = 1
        for tri in inmesh.triangles:
            nodes = [i + 1 for i in tri.nodes] # femmesh format is 1-based, not 0
            print >> outmeshfile, tri_index, tri.material,\
                  nodes[0], nodes[1], nodes[2]
            tri_index += 1

        print >> outmeshfile, "ENDBLOCK"

    def writeTets(self):
        """
        Writes tetraheders to the femmesh output file
        """

        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK tets"
        print >> outmeshfile, len(inmesh.tetraheders)

        # first write out first order tets
        tet_index = 1
        for tet in inmesh.tetraheders:
            nodes = [i + 1 for i in tet.nodes] # femmesh format is 1-based, not 0
            print >> outmeshfile, tet_index, tet.material,\
                  nodes[0], nodes[1], nodes[2], nodes[3]
            tet_index += 1

        print >> outmeshfile, "ENDBLOCK"

    def writeTet2s(self):
        """
        Writes second order tetraheders to the femmesh output file
        """

        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK tet2s"
        print >> outmeshfile, len(inmesh.tetraheders2)

        # now print second order tets
        tet2_index = 1
        for tet2 in inmesh.tetraheders2:
            nodes = [i + 1 for i in tet2.nodes] # femmesh format is 1-based, not 0
            print >> outmeshfile, tet2_index, tet2.material,\
                  nodes[0], nodes[1], nodes[2], nodes[3], nodes[4], nodes[5], nodes[6], nodes[7], nodes[8], nodes[9]
            tet2_index += 1

        print >> outmeshfile, "ENDBLOCK"
    
    def writeNodes(self):
        """
        Writes nodes to the femmesh output file
        """
        outmeshfile = self.outmeshfile
        inmesh = self.inmesh
        print >> outmeshfile, "BLOCK nodes"
        print >> outmeshfile, len(inmesh.nodes)

        node_index = 1
        for node in inmesh.nodes:
            print >> outmeshfile, node_index, node[0], node[1], node[2]
            node_index += 1

        print >> outmeshfile, "ENDBLOCK"
        
    pass

class Mesh(object):

    
    def __init__(self):
        # Ordered lists of x,y,z coords, as [[x1,y1,z1],[x2,y2,z2]] etc.
        self.nodes = []
        self.triangles = []             # Ordered list of Triangle objects
        self.tetraheders = []           # Ordered list of Tetrahedron objects
        self.tetraheders2 = []          # ordered list of second order tetrahedron objects
        self.quadrilaterals = []        # as above, Quadrilaterals
        self.bricks = []
        self.__element_list_map ={
            'triangle': self.triangles,
            'tetrahedron': self.tetraheders,
            'tetrahedron2': self.tetraheders2,
            'quadrilateral': self.quadrilaterals,
            'brick': self.bricks
            }

        pass
    
    def addNode(self, node_coords):
        """
        Adds a node to the end of the node list.

        node_coords should be a list [x y z], where x, y and z are the
        coordinates of the node.
        """
        
        self.nodes.extend([ [float(i) for i in node_coords] ])
        pass
    
    def addElement(self, element):
        """
        Adds an element to the mesh. element should be an instance of
        a class that is derived from Element, and have a shapeid
        attribute known to the class of which this method is a
        member.
        """
        
        self.__element_list_map[element.shapeid].extend([element])
        pass

    def setElementProperties(self, shapeid, element_nos, material):
        """
        Applies property to elements of type shapeid enumerated in
        list element_nos
        """
        elements = self.__element_list_map[shapeid]
        for i in element_nos:
            elements[i].material = material

        pass
    
    pass




# Still have to decide where to put this!


if __name__=='__main__':

    inmeshfile = open(inmeshfilename, 'rU')
    outmeshfile = open(outmeshfilename, 'w')
    mesh = Mesh()                           # internal mesh data representation
    outmesh = OutputMesh(outmeshfile, mesh) # FEMFEKO femmesh output object
    #inmesh = InputGmshMesh(inmeshfile, mesh)# object imesh is of class InputGmshMesh() or InputUNVMesh() depending on the type of file to be converted.
    inmesh = InputGmsh2Mesh(inmeshfile, mesh)
    #inmesh = InputUNVMesh( inmeshfile, mesh)
    #inmesh = InputNetgenMesh(inmeshfile, mesh)
    inmesh.readMesh()
    outmesh.writeMesh()
    outmeshfile.close()
        
