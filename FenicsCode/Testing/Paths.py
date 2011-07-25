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
# Authors:
# Neilen Marais <nmarais@gmail.com>
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

def get_module_path_file(filename, __file__, mode='r'):
    """Return file object located in current module path as specified by __file__
    """
    return open(get_module_path_filename(filename, __file__), mode)
