import copy
import cmath
import random
import numpy as np


def laguerre_method(polynomial, z, tolerance=1e-8):
    degree = polynomial.o
    first_deriv = polynomial.deriv()
    second_deriv = first_deriv.deriv()
    while True:
        P_z = polynomial(z)
        first_deriv_z = first_deriv(z)
        second_deriv_z = second_deriv(z)

        reduced_degree = degree - 1
        sqrt = cmath.sqrt(reduced_degree * (reduced_degree * (first_deriv_z ** 2) - degree * P_z * second_deriv_z))
        addition = first_deriv_z + sqrt
        subtraction = first_deriv_z - sqrt
        denominator = addition if abs(addition) > abs(subtraction) and abs(addition) != 0 else subtraction
        z_new = z - degree * P_z / denominator
        if abs(z_new - z) < tolerance:
            return z_new
        z = z_new


def special_forward_substitution(A, vec_b):
    size = len(vec_b)
    x = []
    x.append(vec_b[0])
    for i in range(1, size):
        x.append(vec_b[i] - A[i, i - 1] * x[i - 1])

    return x


def deflate(polynomial, root):
    vec_b = polynomial.coef[:-1]
    size = len(vec_b)
    A = np.eye(size) - np.eye(size, k=-1) * root
    x = special_forward_substitution(A, vec_b)

    return np.poly1d(x)


def solve_quadratic_equation(quadratic_equation):
    a, b, c = quadratic_equation
    sqrt_delta = cmath.sqrt(b ** 2 - 4 * a * c)
    x1 = (-b + sqrt_delta) / (2 * a)
    x2 = (-b - sqrt_delta) / (2 * a)
    return x1, x2


def solve(polynomial):
    roots = []
    polynomial_copy = copy.deepcopy(polynomial)
    z = random.random()
    while polynomial_copy.o > 2:
        smooth_root = laguerre_method(polynomial, laguerre_method(polynomial_copy, z))
        roots.append(smooth_root)
        polynomial_copy = deflate(polynomial_copy, smooth_root)

    z1, z2 = solve_quadratic_equation(polynomial_copy)

    smooth_z1 = laguerre_method(polynomial, z1)
    roots.append(smooth_z1)

    smooth_z2 = laguerre_method(polynomial, z2)
    roots.append(smooth_z2)
    return roots

def print_solution(polynomial, roots):
    print("Calcualated roots for polynomial:")
    print(polynomial)
    print("are:", roots)
    print()

first_polynomial = np.poly1d([243, -486, 783, -990, 558, -28, -72, 16], variable='z')
second_polynomial = np.poly1d([1, 1, 3, 2, -1, -3, -11, -8, -12, -4, -4], variable='z')
third_polynomial = np.poly1d([1, 1j, -1, -1j, 1], variable='z')

print_solution(first_polynomial, solve(first_polynomial))
print_solution(second_polynomial, solve(second_polynomial))
print_solution(third_polynomial, solve(third_polynomial))