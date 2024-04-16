import numpy as np
import matplotlib.pyplot as plt


def function(x):
    return (x ** 4) / 4 - (x ** 2) / 2 - x / 16


def narrowing_strategy(a, b, c, d):
    fun_a = function(a)
    fun_b = function(b)
    fun_c = function(c)
    fun_d = function(d)

    if fun_d < fun_b:
        if d < b:
            c = b
            fun_c = fun_b
            b = d
            fun_b = fun_d
        else:
            a = b
            fun_a = fun_b
            b = d
            fun_b = fun_d
    else:
        if d < b:
            a = d
            fun_a = fun_d
        else:
            c = d
            fun_c = fun_d

    return a, b, c, fun_a, fun_b, fun_c


def golden_ratio_method(a, c, w=0.381966, max_iterations=100000, tolerance=1e-6):
    b = (a + c) / 2
    distances = []
    iteration = 0

    while True:
        b_prev = b
        abs_b_a = abs(b - a)
        abs_c_b = abs(c - b)
        if abs_b_a > abs_c_b:
            d = a + w * abs_b_a
        else:
            d = b + w * abs_c_b

        a, b, c, fun_a, fun_b, fun_c = narrowing_strategy(a, b, c, d)
        iteration += 1
        abs_c_a = abs(c - a)
        distances.append(abs_c_a)

        if abs_c_a < tolerance * (abs(b_prev) + abs(d)):
            break

        if iteration == max_iterations:
            print("max iterations reached")
            break

    return fun_b, iteration, distances


def calculate_d(a, b, c, fun_a, fun_b, fun_c):
    a_times_fc_fb = a * (fun_c - fun_b)
    b_times_fa_fc = b * (fun_a - fun_c)
    c_times_fb_fa = c * (fun_b - fun_a)

    d = (a_times_fc_fb * a + b_times_fa_fc * b + c_times_fb_fa * c) / (a_times_fc_fb + b_times_fa_fc + c_times_fb_fa)
    d /= 2

    return d


def brent_method(a, c, max_iterations=100000, tolerance=1e-6):
    b = (a + c) / 2
    distances = []
    iteration = 0

    while True:
        d = calculate_d(a, b, c, function(a), function(b), function(c))

        b_prev = b
        abs_a_c_prev = abs(a - c)
        a, b, c, fun_a, fun_b, fun_c = narrowing_strategy(a, b, c, d)

        if a >= d or d >= c or abs(a - c) >= abs_a_c_prev / 2:
            d = (a + c) / 2
            a, b, c, fun_a, fun_b, fun_c = narrowing_strategy(a, b, c, d)

        iteration += 1
        abs_c_a = abs(c - a)
        distances.append(abs_c_a)

        if abs_c_a < tolerance * (abs(b_prev) + abs(d)):
            break

        if iteration == max_iterations:
            print("max iterations reached")
            break


    return fun_b, iteration, distances


first_a = 0
first_c = 2

golden_min, golden_iterations, golden_distances = golden_ratio_method(first_a, first_c)
print("golden ratio method min:", golden_min)

brent_min, brent_iterations, brent_distances = brent_method(first_a, first_c)
print("Brent method min: ", brent_min)

plt.plot(np.arange(0, golden_iterations), golden_distances, linewidth=0.5, color="#C7234F")
plt.scatter(np.arange(0, golden_iterations), golden_distances, label='golden_ratio_method()', s=1, color="#C7234F")

plt.plot(np.arange(0, brent_iterations), brent_distances, linewidth=0.5, color="#4318BA")
plt.scatter(np.arange(0, brent_iterations), brent_distances, label='brent_method()', s=1, color="#4318BA")

plt.xlabel('Iteration')
plt.ylabel('|c-a|')
plt.legend()

plt.savefig('wykres.png', dpi=500)
