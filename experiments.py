import random
import time

from matrices import random_matrix, random_vector, hilbert_matrix, mat_vec
from metrics import residual_norm, relative_error
from solvers import (
    gaussian_no_pivot,
    gaussian_partial_pivot,
    lu_decomposition,
    solve_with_lu,
)
from timing import measure_time


def print_table(headers, rows):
    headers = [str(h) for h in headers]
    rows = [[str(value) for value in row] for row in rows]

    widths = [len(h) for h in headers]
    for row in rows:
        for i in range(len(row)):
            widths[i] = max(widths[i], len(row[i]))

    print(" | ".join(headers[i].ljust(widths[i]) for i in range(len(headers))))
    print("-+-".join("-" * w for w in widths))
    for row in rows:
        print(" | ".join(row[i].ljust(widths[i]) for i in range(len(row))))


# 4.1 - время решения одной системы разными методами
def experiment_single_system(sizes, seed=42):
    print("\n=== 4.1 Сравнение времени решения одной системы ===")
    headers = ["n", "Gauss", "Gauss pivot", "LU decomp", "LU solve",
               "LU total", "Residual (pivot)", "Residual (LU)"]
    rows = []
    rng = random.Random(seed)

    for n in sizes:
        print(f"Размер n = {n}")
        a = random_matrix(n, rng)
        b = random_vector(n, rng)

        try:
            _, t_gauss = measure_time(gaussian_no_pivot, a, b)
            t_gauss = f"{t_gauss:.6f}"
        except ValueError as error:
            t_gauss = f"error: {error}"

        x_pivot, t_pivot = measure_time(gaussian_partial_pivot, a, b)
        (l, u), t_lu_decomp = measure_time(lu_decomposition, a)
        x_lu, t_lu_solve = measure_time(solve_with_lu, l, u, b)

        rows.append([
            n,
            t_gauss,
            f"{t_pivot:.6f}",
            f"{t_lu_decomp:.6f}",
            f"{t_lu_solve:.6f}",
            f"{t_lu_decomp + t_lu_solve:.6f}",
            f"{residual_norm(a, x_pivot, b):.3e}",
            f"{residual_norm(a, x_lu, b):.3e}",
        ])

    print_table(headers, rows)


# 4.2 - выгода LU при нескольких правых частях
def experiment_multiple_rhs(n, k_values, seed=42):
    print("\n=== 4.2 Экономия времени при множественных правых частях ===")
    headers = ["k", "Gauss pivot total", "LU decomp", "LU all solves", "LU total"]
    rows = []
    rng = random.Random(seed)
    a = random_matrix(n, rng)

    (l, u), t_lu_decomp = measure_time(lu_decomposition, a)

    for k in k_values:
        print(f"Количество правых частей k = {k}")
        right_parts = [random_vector(n, rng) for _ in range(k)]

        start = time.perf_counter()
        for b in right_parts:
            gaussian_partial_pivot(a, b)
        t_gauss_total = time.perf_counter() - start

        start = time.perf_counter()
        for b in right_parts:
            solve_with_lu(l, u, b)
        t_lu_solves = time.perf_counter() - start

        rows.append([
            k,
            f"{t_gauss_total:.6f}",
            f"{t_lu_decomp:.6f}",
            f"{t_lu_solves:.6f}",
            f"{t_lu_decomp + t_lu_solves:.6f}",
        ])

    print_table(headers, rows)


# 4.3 - точность на матрицах Гильберта
def experiment_hilbert(sizes):
    print("\n=== 4.3 Проверка точности на матрицах Гильберта ===")
    headers = ["n", "Method", "Relative error", "Residual"]
    rows = []

    methods = [
        ("Gauss", gaussian_no_pivot),
        ("Gauss pivot", gaussian_partial_pivot),
    ]

    for n in sizes:
        h = hilbert_matrix(n)
        x_exact = [1.0] * n
        b = mat_vec(h, x_exact)

        for name, method in methods:
            try:
                x_approx = method(h, b)
                rows.append([n, name,
                             f"{relative_error(x_approx, x_exact):.3e}",
                             f"{residual_norm(h, x_approx, b):.3e}"])
            except ValueError as error:
                rows.append([n, name, f"error: {error}", "-"])

        try:
            l, u = lu_decomposition(h)
            x_lu = solve_with_lu(l, u, b)
            rows.append([n, "LU",
                         f"{relative_error(x_lu, x_exact):.3e}",
                         f"{residual_norm(h, x_lu, b):.3e}"])
        except ValueError as error:
            rows.append([n, "LU", f"error: {error}", "-"])

    print_table(headers, rows)
