import csv, os
#import driver

data_path = 'analytical'
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
    E_data[key] = list([float(x) for x in l]
                       for l in csv.reader(file(filename)))
    
problem_data = {}
for l in file(os.path.join(data_path, 'problem_data.txt')):
    k, v = l.split()
    v = float(v)
    problem_data[k]=v
    
