import copy
import numpy as np

np.set_printoptions(suppress=True)


def householder(A):
    size = len(A)
    R = copy.deepcopy(A)
    for i in range(size - 1):
        alpha = -np.sign(R[i + 1][i]) * np.linalg.norm(R[i + 1:, i])
        r = np.sqrt((alpha * alpha - R[i + 1, i] * alpha) / 2)
        v = np.zeros_like(R[i:, i])
        v[1] = (R[i + 1, i] - alpha) / (2 * r)
        v[2:] = (R[i + 2:, i] / (2 * r))
        P = np.eye(len(v)) - 2 * np.outer(v, v)
        R[i:, i:] = np.dot(np.dot(P, R[i:, i:]), P)
    return R


def my_givens(A, i, j): #Uwaga uzycie modyfikuje maciez A (ze wzgledow optymalizacji szybkosci operacji)
    xi = A[i, i]
    xj = A[j, i]
    r = np.sqrt(xi * xi + xj * xj)
    c = xi / r
    s = xj / r

    size = len(A)
    G = np.eye(size)
    G[i, i] = c
    G[i, j] = s
    G[j, i] = -s
    G[j, j] = c

    A[i: j + 2, i: j + 2] = np.dot(G[i: j + 2, i: j + 2], A[i: j + 2, i: j + 2])
    return A, G


def qr_tridiagonal(A_tridiagonal):
    size = len(A_tridiagonal)
    R = copy.deepcopy(A_tridiagonal)
    QT = np.eye(size)
    for i in range(size - 1):
        R, G = my_givens(R, i, i + 1)
        QT = np.dot(G, QT)
    return QT.T, R


def eigenvalues_qr(A_tridiagonal, iterations=None, tolerance=1e-10):
    size = len(A_tridiagonal)
    if iterations is None:
        iterations = size * size
    Q, R = qr_tridiagonal(A_tridiagonal)
    for i in range(iterations):
        A_new = np.dot(R, Q)
        Q, R = qr_tridiagonal(A_new)
        if np.abs(A_new[np.triu_indices(size, 1)]).max() < tolerance:
            break
    return np.diag(A_new)


A = np.array([[19, 13, 10, 10, 13, -17],
              [13, 13, 10, 10, -11, 13],
              [10, 10, 10, -2, 10, 10],
              [10, 10, -2, 10, 10, 10],
              [13, -11, 10, 10, 13, 13],
              [-17, 13, 10, 10, 13, 19]]) / 12

A_tridiagonal = householder(A)
print("A_tridiagonal:")
print(A_tridiagonal)
print("")
eigenvalues = eigenvalues_qr(A_tridiagonal)
print("Eigenvalues =", eigenvalues)
