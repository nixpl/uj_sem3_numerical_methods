import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)


def function(x):
    return 1 / (1 + 5 * x * x)


def calculate_y_array(x_array):
    y_array = []
    for x in x_array:
        y_array.append(function(x))

    return y_array


def lagrange_interpolation(x, y):
    n = len(x)
    coefficients = np.zeros_like(x, dtype=float)

    for i in range(n):
        numerator = np.poly1d([1])
        denominator = 1

        for j in range(n):
            if j != i:
                numerator *= np.poly1d([1, -x[j]])
                denominator *= (x[i] - x[j])

        coefficients += (numerator / denominator) * y[i]

    return coefficients


x = np.array([-7/8, -5/8, -3/8, -1/8, 1/8, 3/8, 5/8, 7/8])
y = calculate_y_array(x)
polynomial_coefficients = lagrange_interpolation(x, y)
print("wielomian interpolacyjny:")
print(np.poly1d([round(coefficient, 4) for coefficient in polynomial_coefficients]))
# print(np.poly1d(polynomial_coefficients))

x_plt = np.linspace(min(x), max(x))
y_plt = np.polyval(polynomial_coefficients, x_plt)
plt.plot(x_plt, y_plt, color = "GREEN", label="Interpolation polynomial")
plt.scatter(x, y, color = "PURPLE", label="Calculated points")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()
plt.savefig('wykres.png')
