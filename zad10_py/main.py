import numpy as np
import matplotlib.pyplot as plt


def function(x):
    return 1 / (1 + 5 * x * x)


def calculate_y_array(x_array):
    y_array = []
    for x in x_array:
        y_array.append(function(x))

    return y_array


def calculate_weights(size, d=3):
    weights = []
    factorials = {}

    def get_factorial(n):
        if n == 0 or n == 1:
            return 1

        if n in factorials:
            return factorials[n]

        else:
            fac = n * get_factorial(n - 1)
            factorials[n] = fac
            return fac


    fac_d = get_factorial(d)

    for k in range(size + 1):
        sign = 1

        if (k - d) % 2 == 1:
            sign = -1

        summation = 0
        for i in range(max(0, k - d), min(size - d + 1, k + 1)):
            summation += fac_d / (get_factorial(k - i) * get_factorial(d - k + i))

        weight = sign * summation
        weights.append(weight)

    return weights


def fh_interpolation_function(x, weights, x_array, y_array):  # r(x)
    if x in x_array:
        return function(x)

    denominator = 0
    numerator = 0
    iterations = len(x_array)
    for k in range(iterations):
        fraction = weights[k] / (x - x_array[k])
        denominator += fraction
        numerator += fraction * y_array[k]

    r = numerator / denominator

    return r


x_array = np.array([-7 / 8, -5 / 8, -3 / 8, -1 / 8, 1 / 8, 3 / 8, 5 / 8, 7 / 8])
y_array = calculate_y_array(x_array)

weights = calculate_weights(len(x_array))

plt.scatter(x_array, y_array, color="BLACK", label="Calculated points")

x_plt = np.linspace(min(x_array), max(x_array), 1000)
y_plt = [fh_interpolation_function(x, weights, x_array, y_array) for x in x_plt]

plt.plot(x_plt, y_plt, label=f"Interpolation")
plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.legend(loc='lower center', fontsize='x-small')

plt.savefig('wykres.png')
