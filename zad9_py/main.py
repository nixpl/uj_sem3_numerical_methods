import numpy as np
import math
import matplotlib.pyplot as plt

np.set_printoptions(suppress=True)


def function(x):
    return 1 / (1 + 5 * x * x)


def calculate_y_array(x_array):
    y_array = []
    for x in x_array:
        y_array.append(function(x))

    return y_array


def create_matrix_A_for_equal_distance(size):
    size -= 2
    A = 4 * np.eye(size) + np.eye(size, k=-1) + np.eye(size, k=1)
    return A


def create_vector_b_for_equal_distance(y, distance):
    iterations = np.size(y) - 2
    multiplier = 6 / (distance ** 2)
    b = []
    for i in range(iterations):
        b.append(multiplier * (y[i] - 2 * y[i + 1] + y[i + 2]))

    return b


def cholseky_tridiag(A):
    C = np.zeros_like(A)
    iterations = len(A)

    C[0, 0] = math.sqrt(A[0, 0])
    for i in range(iterations):
        buff = A[i, i]
        for j in range(i):
            buff -= C[i, j] ** 2

        C[i, i] = math.sqrt(buff)

        if i + 1 == iterations:
            break

        C[i + 1, i] = A[i + 1, i] / C[i, i]

    return C


def solve_for_cholesky(C, b):
    def solve_C(C, b):
        iterations = np.size(b)
        x = []
        x.append(b[0] / C[0, 0])
        for i in range(1, iterations):
            x.append((b[i] - C[i, i - 1] * x[i - 1]) / C[i][i])

        return x

    def solve_CT(C, b):
        iterations = np.size(b) - 1
        x = np.zeros_like(b)
        x[iterations] = b[iterations] / C[iterations, iterations]
        iterations -= 1
        for i in range(iterations, -1, -1):
            x[i] = (b[i] - C[i + 1, i] * x[i + 1]) / C[i][i]

        return x

    return solve_CT(C, solve_C(C, b))


def calculate_cubic_coefficients(x):
    A = create_matrix_A_for_equal_distance(np.size(x))
    b = create_vector_b_for_equal_distance(y, x[1] - x[0])
    C = cholseky_tridiag(A)
    xis = []
    xis.append(0)
    xis.extend(solve_for_cholesky(C, b))
    xis.append(0)

    return xis


def calculate_cubic_function(section_num, x_array, y_array, xis):
    distance = x_array[1] - x_array[0]
    a = np.poly1d([-1, x_array[section_num + 1]]) / distance
    b = np.poly1d([1, -x_array[section_num]]) / distance
    c = (((a ** 3) - a) * distance ** 2) / 6
    d = (((b ** 3) - b) * distance ** 2) / 6

    polynomial = a * y_array[section_num] + b * y_array[section_num + 1] + c * xis[section_num] + d * \
                 xis[section_num + 1]
    return polynomial


def create_cubic_plot(x, y, xis):
    plot_sections_count = np.size(x) - 1
    for i in range(plot_sections_count):
        x_plt = np.linspace(x[i], x[i + 1])
        y_plt = np.polyval(calculate_cubic_function(i, x, y, xis), x_plt)
        color = plt.cm.rainbow(i / plot_sections_count)
        plt.plot(x_plt, y_plt, color=color, label=f"Interpolation for [{x[i]}, {x[i + 1]}]")


x = np.array([-7 / 8, -5 / 8, -3 / 8, -1 / 8, 1 / 8, 3 / 8, 5 / 8, 7 / 8])
y = calculate_y_array(x)
xis = calculate_cubic_coefficients(x)
print("wspolczynniki kubiczne:", xis)

plt.scatter(x, y, color = "BLACK", label="Calculated points")
create_cubic_plot(x, y, xis)

plt.xlabel("x")
plt.ylabel("y")
plt.grid()
plt.legend(loc='lower center', fontsize='x-small')

plt.savefig('wykres.png')
