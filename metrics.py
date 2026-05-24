from matrices import mat_vec, vector_sub, norm2


# невязка - проверяет, хорошо ли выполняется уравнение Ax = b
def residual_norm(a, x, b):
    return norm2(vector_sub(mat_vec(a, x), b))


# относительная погрешность - насколько найденное решение близко к точному
def relative_error(x_approx, x_exact):
    numerator = norm2(vector_sub(x_approx, x_exact))
    denominator = norm2(x_exact)
    if denominator == 0:
        return numerator
    return numerator / denominator
