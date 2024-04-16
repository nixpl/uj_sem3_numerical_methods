import numpy as np
np.set_printoptions(suppress=True, linewidth=200)

def qr_shift_algorithm(A, max_iterations=10000, tolerance=1e-10, shift=1e-1):
    size = len(A)
    eigen_vectors = np.eye(size)
    shifted_eye = shift * np.eye(size)
    for i in range(max_iterations):
        Q, R = np.linalg.qr(A - shifted_eye)
        A = np.dot(R, Q) + shifted_eye
        eigen_vectors = np.dot(eigen_vectors, Q)

        if np.abs(A[np.triu_indices(size, 1)]).max() < tolerance:
            break

    eigen_values = np.diag(A)

    return eigen_values, eigen_vectors

def to_complex_vectors(noncomplex_vectors):
    half_size = len(noncomplex_vectors) // 2
    re = noncomplex_vectors[:half_size, :]
    im = noncomplex_vectors[half_size:, :]
    complex_vectors = re + im * 1j

    return complex_vectors


A = np.array([[0, 1, 0, 0],
              [1, 0, 0, 0],
              [0, 0, 0, 1],
              [0, 0, 1, 0]])

B = np.array([[0, 0, 0, -1],
              [0, 0, -1, 0],
              [0, 1, 0, 0],
              [1, 0, 0, 0]])

H = np.block([[A, -B], [B, A]])

eigen_values, eigen_vectors = qr_shift_algorithm(H)


unique_eigen_values, unique_indices = np.unique(eigen_values, return_index=True)
unique_eigen_vectors = eigen_vectors[:, unique_indices]

unique_complex_eigen_vectors = to_complex_vectors(unique_eigen_vectors)

print(H)
print("")

print(unique_eigen_values)
print("")

print(unique_complex_eigen_vectors)
print("")
