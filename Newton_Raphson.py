def find_distance_newton(x0, y0, f, df, ddf, initial_guess=0.0,
tolerance=1e-7, max_iter=100):
    x = initial_guess
    for _ in range(max_iter):
        # D'(x)
        D_prime = 2 * (x - x0) + 2 * (f(x) - y0) * df(x)
        # D''(x)
        D_double_prime = 2 + 2 * (df(x)**2) + 2 * (f(x) - y0) * ddf(x)
# Newton-Raphson update step
        next_x = x - D_prime / D_double_prime
        if abs(next_x - x) < tolerance:
            break
        x = next_x
        shortest_distance = ((x - x0)**2 + (f(x) - y0)**2)**0.5
    return shortest_distance, x