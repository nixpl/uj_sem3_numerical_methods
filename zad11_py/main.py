import numpy as np
import math
np.set_printoptions(suppress=True)

def function(x):
    return math.sin(math.pi * (1 + math.sqrt(x)) / (1 + x ** 2)) * math.exp(-x)


def calculate_b(x, tolerance=1e-7):
    while math.exp(-x) > tolerance:
        x += 1
    return x


def trapeze_method(points, sum, N):
    end = len(points) - 1
    new_points = []
    for i in range(end):
        p1 = points[i]
        p3 = points[i+1]
        p2 = (p1 + p3)/2
        new_points.append(p1)
        new_points.append(p2)

    new_points.append(points[-1])

    end = len(new_points) - 1
    for i in range(1, end, 2):
        sum += function(new_points[i])

    h = new_points[1] - new_points[0]
    I = h * sum
    return I, sum, new_points


def calculate_integral_romberg(a, b, tolerance=1e-7):

    A_last_row = []
    A_current_row = []

    k = 0
    N = 2 ** k
    points = [a, b]
    sum = (function(a) + function(b)) / 2
    I, sum, points = trapeze_method(points, sum, N)
    A_last_row.append(I)

    while True:
        k += 1
        N = 2 ** k
        I, sum, points = trapeze_method(points, sum, N)
        A_current_row.append(I)

        end = len(A_last_row) + 1
        for i in range(1, end):
            power = 4 ** i
            I = 1 / (power - 1) * (power * A_current_row[i-1] - A_last_row[i-1])
            A_current_row.append(I)

        if abs(A_current_row[-1] - A_last_row[-1]) < tolerance:
            break

        A_last_row = A_current_row.copy()
        A_current_row.clear()

    return I


a = 0
b = calculate_b(a)

I = calculate_integral_romberg(a, b)

print("Integral:", I)
