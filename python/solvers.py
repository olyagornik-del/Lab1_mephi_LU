from matrices import copy_matrix, copy_vector

EPS = 1e-12


# решаем Ux = y, идём снизу вверх
def backward_substitution(u, y):
    n = len(u)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(u[i][j] * x[j] for j in range(i + 1, n))
        pivot = u[i][i]
        if abs(pivot) < EPS:
            raise ValueError(f"Нулевой или слишком маленький диагональный элемент U[{i}][{i}]")
        x[i] = (y[i] - s) / pivot
    return x


# решаем Ly = b, идём сверху вниз
def forward_substitution(l, b):
    n = len(l)
    y = [0.0] * n
    for i in range(n):
        s = sum(l[i][j] * y[j] for j in range(i))
        pivot = l[i][i]
        if abs(pivot) < EPS:
            raise ValueError(f"Нулевой или слишком маленький диагональный элемент L[{i}][{i}]")
        y[i] = (b[i] - s) / pivot
    return y


# классический метод Гаусса (без выбора ведущего элемента)
def gaussian_no_pivot(a, b):
    a = copy_matrix(a)
    b = copy_vector(b)
    n = len(a)

    for k in range(n - 1):
        pivot = a[k][k]
        if abs(pivot) < EPS:
            raise ValueError(f"Нулевой или слишком маленький ведущий элемент на шаге {k}")
        for i in range(k + 1, n):
            factor = a[i][k] / pivot
            a[i][k] = 0.0
            for j in range(k + 1, n):  # пересчитываем остальные элементы строки
                a[i][j] -= factor * a[k][j]
            b[i] -= factor * b[k]

    return backward_substitution(a, b)


# метод Гаусса с частичным выбором ведущего элемента по столбцу
def gaussian_partial_pivot(a, b):
    a = copy_matrix(a)
    b = copy_vector(b)
    n = len(a)

    for k in range(n - 1):
        pivot_row = max(range(k, n), key=lambda i: abs(a[i][k]))  # ищем макс эл в столбце
        if abs(a[pivot_row][k]) < EPS:
            raise ValueError(f"Столбец {k} вырожден: ведущий элемент слишком мал")
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            b[k], b[pivot_row] = b[pivot_row], b[k]

        pivot = a[k][k]
        for i in range(k + 1, n):
            factor = a[i][k] / pivot
            a[i][k] = 0.0
            for j in range(k + 1, n):
                a[i][j] -= factor * a[k][j]
            b[i] -= factor * b[k]

    return backward_substitution(a, b)


# LU-разложение без перестановок: A = L * U
def lu_decomposition(a):
    n = len(a)
    l = [[0.0] * n for _ in range(n)]
    u = [[0.0] * n for _ in range(n)]
    for i in range(n):
        l[i][i] = 1.0

    for k in range(n):
        for j in range(k, n):
            u[k][j] = a[k][j] - sum(l[k][s] * u[s][j] for s in range(k))
            # тк a_kj = l_k0 * u_0j + l_k1 * u_1j + ... + l_kk-1 * u_k-1j + l_kk * u_kj
            # a_kj знаем, l_kk равен 1, те чтобы найти u_kj, из a_kj вычитаем сумму

        if abs(u[k][k]) < EPS:
            raise ValueError(f"LU-разложение без перестановок невозможно: маленький pivot на шаге {k}")

        for i in range(k + 1, n):
            l[i][k] = (a[i][k] - sum(l[i][s] * u[s][k] for s in range(k))) / u[k][k]

    return l, u


def solve_with_lu(l, u, b):
    y = forward_substitution(l, b)
    return backward_substitution(u, y)
