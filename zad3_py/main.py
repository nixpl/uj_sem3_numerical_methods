import copy
import numpy as np
np.set_printoptions(suppress=True)

def find_eigens(A, quantity, iterations=1000, tolerance=1e-10):
    size = len(A)
    eigons = []

    for k in range(quantity):
        vec = np.random.rand(size)
        vec /= np.linalg.norm(vec)
        for j in range(iterations):
            Avec = np.dot(A, vec)
            vec_new = Avec / np.linalg.norm(Avec)

            if np.linalg.norm(vec_new - vec) < tolerance:
                break

            vec = vec_new

        eigenval = np.dot(vec, np.dot(A, vec)) / np.dot(vec, vec)
        eigon = (eigenval, vec)
        eigons.append(eigon)
        A -= eigenval * np.outer(vec, vec)

    return eigons


def print_egons(eigens):
    for eigen in eigens:
        print("Wartosc wlasna:", eigen[0])
        print("Odpowiadajacy jej wektor wlasny:", eigen[1])
        print("")

A = np.array([[19, 13, 10, 10, 13, -17],
              [13, 13, 10, 10, -11, 13],
              [10, 10, 10, -2, 10, 10],
              [10, 10, -2, 10, 10, 10],
              [13, -11, 10, 10, 13, 13],
              [-17, 13, 10, 10, 13, 19]]) / 12

# tworze kopie na potrzeby testow
A_copy = copy.deepcopy(A)

eigens = find_eigens(A, 2)
print_egons(eigens)

