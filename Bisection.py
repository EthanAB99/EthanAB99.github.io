import math


def golden_section_search(x0, y0, f, a, b, tolerance=1e-7):
    # a and b should be chosen inside the domain of f
    phi = (1 + math.sqrt(5)) / 2
    resphi = 2 - phi
    x1 = a + resphi * (b - a)
    x2 = b - resphi * (b - a)

    def dist_sq(x):
        return (x - x0) ** 2 + (f(x) - y0) ** 2

    f_x1 = dist_sq(x1)
    f_x2 = dist_sq(x2)

    while abs(b - a) > tolerance:
        if f_x1 < f_x2:
            b = x2
            x2 = x1
            f_x2 = f_x1
            x1 = a + resphi * (b - a)
            f_x1 = dist_sq(x1)
        else:
            a = x1
            x1 = x2
            f_x1 = f_x2
            x2 = b - resphi * (b - a)
            f_x2 = dist_sq(x2)

    best_x = (a + b) / 2
    return math.sqrt(dist_sq(best_x)), best_x
