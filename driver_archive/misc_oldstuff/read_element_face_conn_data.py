"""
Parses eMAGUS .out file for element face connectivity data
"""

import re

def get_block(file, regex_beg):
    regex_end = re.compile(r'\s+$')
    for line in file:
        if regex_beg.search(line):
            print file.next()
            for line in file:
                if regex_end.match(line):
                    raise StopIteration
                yield line
    pass
outfile = 'NewCode/tests/testdata/test_eMAGUSImport/test_tet.out'
elblock_re = re.compile("ELEMENT FACE INTERCONNECTIVITY DATA")
faceblock_re = re.compile("LIST OF FACES")
edgeblock_re = re.compile("LIST OF EDGES")
numbers_re = re.compile('\d+')


edges = {}
edge_re = re.compile(r'(\d+)(\s+\d+){2}\s+([TF])')
for edge in get_block(file(outfile, 'r'), edgeblock_re):
    m = edge_re.search(edge)
    if m: edges[int(m.group(1))] = (m.group(3) == 'T')






faces = {}
face_re = re.compile(r'(\d+)(\s+\d){3,}\s+([TF])')
for face in get_block(file(outfile, 'r'), faceblock_re):
    m = face_re.search(face)
    if m: faces[int(m.group(1))] = (m.group(3) == 'T')



els = {}                                # Use a dict to kill the element no dupes
for el in get_block(file(outfile, 'r'), elblock_re):
    nos = [int(x) for x in numbers_re.findall(el)]
    elno = nos[0]
    els[elno] = {'connect2elem' : nos[2:6],
                 'connect2face' : nos[6:11]}
    
els = [x[1] for x in sorted(els.iteritems())]
connect2elem = [el['connect2elem'] for el in els]
connect2face = [el['connect2face'] for el in els]
