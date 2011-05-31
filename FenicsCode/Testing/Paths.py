import os

def get_module_path(modfile):
    """Return the current module path, given __file__

    To find out in which directory on the filesystem the currently
    executing code is located, the __file__ special is passed to
    get_module_path(). The path of the current module is returned."""

    return os.path.dirname(modfile)

def get_module_path_filename(filename, __file__):
    """Return filename located in current module path as specified by __file__
    """
    dirname = get_module_path(__file__)
    return os.path.join(dirname, filename)

def get_module_path_file(filename, __file__, mode='f'):
    """Return file object located in current module path as specified by __file__
    """
    return open(get_module_path_filename(filename, __file__), mode)
