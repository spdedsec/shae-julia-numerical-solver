# SHAE Julia Numerical Solver
# Authentic Julia syntax - numerical methods and root finding

"""
    bisection(f, a, b, tol=1e-6, max_iter=100)

Find root of f(x) = 0 using bisection method on [a, b].
"""
function bisection(f, a, b, tol=1e-6, max_iter=100)
    if f(a) * f(b) > 0
        error("f(a) and f(b) must have opposite signs")
    end
    for i in 1:max_iter
        c = (a + b) / 2
        if abs(f(c)) < tol || (b - a) / 2 < tol
            return c
        end
        if f(a) * f(c) < 0
            b = c
        else
            a = c
        end
    end
    return (a + b) / 2
end

"""
    newton_raphson(f, f_prime, x0, tol=1e-6, max_iter=100)

Find root using Newton-Raphson method.
"""
function newton_raphson(f, f_prime, x0, tol=1e-6, max_iter=100)
    x = x0
    for i in 1:max_iter
        fx = f(x)
        fpx = f_prime(x)
        if abs(fpx) < 1e-12
            error("Derivative near zero")
        end
        x_new = x - fx / fpx
        if abs(x_new - x) < tol
            return x_new
        end
        x = x_new
    end
    return x
end

"""
    integrate(f, a, b, n=1000)

Numerical integration using trapezoidal rule.
"""
function integrate(f, a, b, n=1000)
    h = (b - a) / n
    s = 0.5 * (f(a) + f(b))
    for i in 1:n-1
        s += f(a + i * h)
    end
    return s * h
end

"""
    solve_linear_system(A, b)

Solve Ax = b using Gaussian elimination.
"""
function solve_linear_system(A, b)
    n = length(b)
    M = [A[i, :] b[i] for i in 1:n]
    for col in 1:n
        # Partial pivoting
        max_row = argmax(abs.(M[col:n, col])) + col - 1
        if M[max_row, col] ≈ 0
            error("Singular matrix")
        end
        M[col, :], M[max_row, :] = M[max_row, :], M[col, :]
        for row in col+1:n
            factor = M[row, col] / M[col, col]
            M[row, :] -= factor * M[col, :]
        end
    end
    x = zeros(n)
    for i in n:-1:1
        x[i] = (M[i, end] - sum(M[i, j] * x[j] for j in i+1:n)) / M[i, i]
    end
    return x
end

function main()
    println("SHAE Julia Numerical Solver")
    println("Available methods:")
    println("  - bisection")
    println("  - newton_raphson")
    println("  - integrate (trapezoidal)")
    println("  - solve_linear_system (Gaussian elimination)")
    println()

    # Example: find root of x^2 - 4 = 0
    f(x) = x^2 - 4
    root = bisection(f, 0.0, 3.0)
    println("Bisection root of x^2 - 4: $root")

    # Example: integrate sin(x) from 0 to pi
    area = integrate(sin, 0.0, π)
    println("Integral of sin(x) from 0 to π: $area")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
