import numpy as np
np.set_printoptions(suppress=True)

def find_eigen_vector(A, eigen_value, iterations=1000, tolerance=1e-10):
    size = len(A)

    vec = np.random.rand(size)
    vec /= np.linalg.norm(vec)
    for j in range(iterations):
        vec_z = np.linalg.solve(A - eigen_value * np.eye(size), vec)
        vec_new = vec_z / np.linalg.norm(vec_z)

        if np.linalg.norm(vec_new - vec) < tolerance:
            break

        vec = vec_new

    return vec_new


A = np.array([[2, -1, 0, 0, 1],
                    [-1, 2, 1, 0, 0],
                    [0, 1, 1, 1, 0],
                    [0, 0, 1, 2, -1],
                    [1, 0, 0, -1, 2]], dtype=np.float64)

EIGEN_VALUE = 0.38197

eigen_vector = find_eigen_vector(A, EIGEN_VALUE)

print("Wektor wlasny dla wartosci wlasnej", EIGEN_VALUE, "to:")
print(eigen_vector)
