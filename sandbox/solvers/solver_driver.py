__author__ = "Evan Lezar"
__date__ = "20 May 2011"

"""
A simple driver to test various LA solvers.
"""

def generate_matrices ( mesh_file ):
    return None, None

def save_matrices ( id, A, b ):
    pass

def load_saved_matrices ():
    return None, None


def main ():
    force_generate = True;
    
    mesh_file = '../../workspace/sphere-r1m-4.femmesh'
    problem_id = 'sphere-r1m-4'
    
    if force_generate:
        A, b = generate_matrices ( mesh_file )
        save_matrices ( problem_id, A, b )
    else:
        A, b = load_matrices ( A, b )
    
    print A
    print b

if __name__ == "__main__":
    main()