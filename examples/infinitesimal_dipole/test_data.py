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
## along with SUCEM. If not, see <http://www.gnu.org/licenses/>. 
##
## Contact: cemagga@gmail.com 
# Authors:
# Neilen Marais <nmarais@gmail.com>
import csv, os
import numpy as N

data_path = 'analytical'
data_path = os.path.join(os.path.dirname(__file__), data_path)
E_data_csv_files = ('E_vals_1_re.csv',
                    'E_vals_1_im.csv',
                    'E_vals_2_re.csv',
                    'E_vals_2_im.csv',
                    'test_pts_1.csv',
                    'test_pts_2.csv',)
E_data = {}
for df in E_data_csv_files:
    key = os.path.splitext(df)[0]
    filename = os.path.join(data_path, df)
    E_data[key] = N.array([[float(x) for x in l]
                       for l in csv.reader(file(filename))], N.float64)
    
problem_data = {}
for l in file(os.path.join(data_path, 'problem_data.txt')):
    k, v = l.split()
    v = float(v)
    problem_data[k]=v
    
E_1 = E_data['E_vals_1_re'] + 1j*E_data['E_vals_1_im']
E_2 = E_data['E_vals_2_re'] + 1j*E_data['E_vals_2_im']
r_1 = E_data['test_pts_1']
r_2 = E_data['test_pts_2']
