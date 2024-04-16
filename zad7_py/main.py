import numpy as np
import matplotlib.pyplot as plt
np.set_printoptions(suppress=True)

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

    coefficients = np.round(coefficients, decimals=4)
    return coefficients


# [n][0] wartosc argumentu x
# [n][1] wartosc funckcji dla argumentu [n][0]
xy = np.array([[-0.75, 1.1309204101562500],
               [-0.50, 2.3203125000000000],
               [-0.25, 1.9284057617187500],
               [0.00, 1.0000000000000000],
               [0.25, 0.0554809570312500],
               [0.50, -0.6015625000000000],
               [0.75, -0.7525024414062500],
               [1.00, 0.0000000000000000]])

coefficients = lagrange_interpolation(xy[:, 0], xy[:, 1])
print("Wspolczynniki w kolejnosci dla zapisywania wielomianu:", coefficients)
x = np.linspace(-1.25, 1.25)
y = np.polyval(coefficients, x)
plt.scatter(xy[:, 0], xy[:, 1], color = "PURPLE", label="Given points")
plt.plot(x, y, color = "GREEN", label="Interpolation polynomial")
plt.xlabel("x")
plt.ylabel("y")
plt.legend()
plt.grid()

plt.savefig('wykres.png')
