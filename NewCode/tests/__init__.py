import TestMeshes

def xfail(testfun):
    def xfail_testfun(*names, **kwargs):
        try:
            testfun(*names, **kwargs)
        except:
            print 'X'
    return xfail_testfun

