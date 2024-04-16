import matplotlib.pyplot as plt
import math

def function(x):
    return math.cos((1 + x) / (x**2 + 0.04)) * math.exp(-(x**2))

def calculate_b(x=0, tolerance=1e-8):
    while math.exp(-(x ** 2)) > tolerance:
        x += 1
    return x

def find_middle_point(a, b):
    return (a+b)/2

def trapeze_method(a, b):
    halfH = (b-a)/2
    return halfH * (function(a) + function(b))

def adaptive_quadrature(a, b, tolerance=1e-8, stack_size = 100):
    stack_overflow_info = False
    mid = find_middle_point(a, b)
    prev_result = trapeze_method(a, b)
    integral = 0
    stack = []

    xplt = []
    yplt = []

    while True:
        left_result = trapeze_method(a, mid)
        right_result = trapeze_method(mid, b)

        if abs(prev_result - (left_result + right_result)) < tolerance:
            integral += prev_result
            xplt.append(mid)
            yplt.append(integral)
            if not stack:
                break

            a, b, prev_result = stack.pop()
            mid = find_middle_point(a, b)

        else:
            if len(stack) == stack_size:
                stack_overflow_info = True
                break

            stack.append([mid, b, right_result])
            # a does not change
            b = mid
            mid = find_middle_point(a, b)
            prev_result = left_result

    return stack_overflow_info, integral, xplt, yplt


b = calculate_b()
a = -b

print("Integral range: ", "[", a, ",", b, "]")

stack_overflow_info, integral, xplt, yplt = adaptive_quadrature(a, b)

if not stack_overflow_info:
    print("Integral:", integral)

    for element in xplt:
        plt.axvline(element, color='grey', linewidth = 0.05)

    plt.plot(xplt, yplt, label="F(x)")

    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend()
    plt.grid()

    plt.savefig("wykres.png", dpi = 1000)
else:
    print("Stack Overflow")




