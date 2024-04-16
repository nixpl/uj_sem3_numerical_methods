import numpy as np
import matplotlib.pyplot as plt


def gauss_seidel(gauss_norms, A, b, x=None, tolerance=1e-8):
    matrix_size = len(b)
    if x == None:
        x = np.zeros(matrix_size)
    norm = tolerance + 1
    while norm > tolerance:
        x_prev = x.copy()
        for i in range(matrix_size):
            sum = 0
            if i >= 4:
                sum += A[i, i - 4] * x[i - 4]
            if i >= 1:
                sum += A[i][i - 1] * x[i - 1]
            if i < matrix_size - 4:
                sum += A[i][i + 4] * x[i + 4]
            if i < matrix_size - 1:
                sum += A[i][i + 1] * x[i + 1]

            x[i] = (e[i] - sum) / A[i, i]

        norm = np.linalg.norm(x - x_prev)
        gauss_norms.append(norm)

    return x


def special_dot(A, p):
    matrix_size = len(A)
    vector = np.zeros(matrix_size)
    for i in range(matrix_size):
        if i >= 4:
            vector[i] += A[i, i - 4] * p[i - 4]
        if i >= 1:
            vector[i] += A[i][i - 1] * p[i - 1]
        if i >= 0:
            vector[i] += A[i][i] * p[i]
        if i < matrix_size - 4:
            vector[i] += A[i][i + 4] * p[i + 4]
        if i < matrix_size - 1:
            vector[i] += A[i][i + 1] * p[i + 1]

    return vector


def conjugate_gradient(gradient_norms, A, b, x=None, tol=1e-8):
    if x is None:
        x = np.zeros(len(b))
    r = b - np.dot(A, x)
    p = r
    matrix_size = len(b)
    i = 0
    stop = False
    while i < matrix_size and stop == False:
        dot_ap = special_dot(A, p)
        dot_rr = np.dot(r, r)
        alpha = dot_rr / np.dot(p, dot_ap)
        r_new = r - alpha * dot_ap
        beta = np.dot(r_new, r_new) / dot_rr
        p_new = r_new + beta * p
        x = x + alpha * p
        p = p_new
        r = r_new

        norm = np.linalg.norm(r_new)
        gradient_norms.append(norm)
        if norm < tol:
            stop = True

    return x


matrixSize = 128
e = np.ones(matrixSize)
A = np.zeros((matrixSize, matrixSize))
for i in range(matrixSize):  # inicjalizacja macierzy A
    A[i, i] = 4
    if i < matrixSize - 1:
        A[i, i + 1] = 1
    if i < matrixSize - 4:
        A[i, i + 4] = 1
    if i >= 1:
        A[i, i - 1] = 1
    if i >= 4:
        A[i, i - 4] = 1

gauss_norms = []
gradient_norms = []

print("Rozwiazanie przy pomocy funcji gauss_seidel(): x =", gauss_seidel(gauss_norms, A, e))

print("Rozwiazanie przy pomocy funcji conjugate_gradient(): x =", conjugate_gradient(gradient_norms, A, e))

x_gauss_values = np.arange(len(gauss_norms))
x_gradient_values = np.arange(len(gradient_norms))

plt.plot(x_gauss_values, gauss_norms, linewidth=0.5, color="#C7234F")
plt.scatter(x_gauss_values, gauss_norms, label='Gauss-Seidel', s=1, color="#C7234F")

plt.plot(x_gradient_values, gradient_norms, linewidth=0.5, color="#4318BA")
plt.scatter(x_gradient_values, gradient_norms, label='Conjugate Gradient', s=1, color="#4318BA")

plt.xlabel('Iteration')
plt.ylabel('Norm')
plt.legend()

plt.savefig('wykres.png')
