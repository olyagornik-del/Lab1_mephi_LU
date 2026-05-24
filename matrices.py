import math


def copy_matrix(a):
    return [list(row) for row in a]


def copy_vector(b):
    return list(b)


def random_matrix(n, rng, low=-1.0, high=1.0):
    return [[rng.uniform(low, high) for _ in range(n)] for _ in range(n)]


def random_vector(n, rng, low=-1.0, high=1.0):
    return [rng.uniform(low, high) for _ in range(n)]


# Hij = 1 / (i + j - 1), если нумеровать с единицы; у нас индексы с нуля
def hilbert_matrix(n):
    return [[1.0 / (i + j + 1) for j in range(n)] for i in range(n)]


# умножение матрицы на вектор
def mat_vec(a, x):
    return [sum(a[i][j] * x[j] for j in range(len(x))) for i in range(len(a))]


def vector_sub(a, b):
    return [a[i] - b[i] for i in range(len(a))]


# длинна вектора
def norm2(v):
    return math.sqrt(sum(value * value for value in v))
