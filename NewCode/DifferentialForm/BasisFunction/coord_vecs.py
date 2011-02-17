import numpy as N

# Reference coordinate-vectors for covariant components
cov_coord_vecs = N.array([[1,0,0,0],
                          [0,1,0,0],
                          [0,0,1,0],
                          [0,0,0,1]], N.float64)

# Reference coordinate-vector "matrix" for contravariant components.  All the
# possible cross products of the covariant component coordinate-vectors are
# represented in a matrix. conv_coord_vecs[i,j] represents grad(lambda_i) X
# grad(lambda_j). Excluding negatives, there are six unique coordinate vector
# possiblities, represented by "unit vectors" of dimension 6.
eye = N.eye(6, dtype=N.float64)
conv_coord_vecs = {}
k = 0
for i in range(3):
    for j in range(i+1,4,1):
        conv_coord_vecs[i,j] = eye[k]
        conv_coord_vecs[j,i] = -eye[k]
        k +=1
