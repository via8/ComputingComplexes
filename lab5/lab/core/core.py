import numpy

from kaucherpy.kaucher import Kaucher as Interval


def immersion(vec, inv=False):
    vec_dim = len(vec) // 2 if inv else len(vec)
    vec_new = numpy.empty(vec_dim if inv else len(vec) * 2, dtype=(Interval if inv else float))
    for j in range(vec_dim):
        if inv:
            vec_new[j] = Interval(-vec[j], vec[j + vec_dim])
        else:
            vec_new[j] = -vec[j].lower
            vec_new[j + vec_dim] = vec[j].upper
    return vec_new


def sign_submatrix(matrix, pos=True):
    result = numpy.zeros(shape=matrix.shape)
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            result[i][j] = max(0.0, matrix[i][j]) if pos else min(0.0, matrix[i][j])
    return result


def calc_initial(A, b):
    M = len(A)
    N = len(A[0])
    A_mid = numpy.array([
        [
            (A[i][j].lower + A[i][j].upper) / 2.0
            for j in range(N)
        ] for i in range(M)
    ])
    blocks_upper = numpy.concatenate((
        sign_submatrix(A_mid), sign_submatrix(A_mid, pos=False)), axis=1)
    blocks_lower = numpy.concatenate((
        sign_submatrix(A_mid, pos=False), sign_submatrix(A_mid)), axis=1)
    imm_matrix = numpy.concatenate((blocks_upper, blocks_lower), axis=0)
    imm_vector = immersion(b)
    return numpy.linalg.solve(imm_matrix, imm_vector)


def subdiff_newton(A, b, initial, epsilon=1e-6):
    def G(x):
        tmp1 = immersion(x, inv=True)
        return immersion(numpy.matmul(A, immersion(x, inv=True))) - immersion(b)

    result = initial
    subdiff = G(result)
    while numpy.linalg.norm(subdiff, 2) > epsilon:
        result -= numpy.matmul(numpy.linalg.inv(subdiff_value(A, result)), subdiff)
        subdiff = G(result)
    return immersion(result, inv=True)


def subdiff_value(A, x):
    dim = len(x)
    value = numpy.zeros(shape=(dim, dim))
    for i in range(dim):
        value[i] = numpy.transpose(subgrad(A, x, i))
    return value


def subgrad(A, x, i):
    dim = len(x) // 2
    value = numpy.zeros(shape=x.shape)
    if 0 <= i < dim:
        for j in range(dim):
            dmax = minor(A, x, i, j)
            tmp1 = max(0.0, A[i][j].lower) * dx(x[j], pos=False) + \
                   min(0.0, A[i][j].upper) * dx(x[j + dim], pos=False) - dmax[0]
            tmp2 = max(0.0, A[i][j].lower) * dx(x[j], pos=False) + \
                   min(0.0, A[i][j].upper) * dx(x[j + dim], pos=False) - dmax[1]
            value[j] -= tmp1
            value[j + dim] -= tmp2
    else:
        for j in range(dim):
            dmax = major(A, x, i - dim, j)
            tmp1 = dmax[0] - \
                   max(0.0, A[i - dim][j].lower) * dx(x[j + dim], pos=False) - \
                   min(0.0, A[i - dim][j].upper) * dx(x[j], pos=False)
            tmp2 = dmax[1] - \
                   max(0.0, A[i - dim][j].lower) * dx(x[j + dim], pos=False) - \
                   min(0.0, A[i - dim][j].upper) * dx(x[j], pos=False)
            value[j] += tmp1
            value[j + dim] += tmp2
    return value


def minor(A, x, i, j):
    dim = len(x) // 2
    tmp1 = max(0.0, A[i][j].upper) * max(0.0, x[j])
    tmp2 = min(0.0, A[i][j].lower) * max(0.0, x[j + dim])
    if tmp1 > tmp2:
        return max(0.0, A[i][j].upper), 0.0
    elif tmp1 == tmp2:
        return 0.5 * max(0.0, A[i][j].upper), 0.5 * min(0.0, A[i][j].lower)
    else:
        return 0.0, min(0.0, A[i][j].lower)


def major(A, x, i, j):
    dim = len(x) // 2
    tmp1 = max(0.0, A[i][j].upper) * max(0.0, x[j + dim])
    tmp2 = min(0.0, A[i][j].lower) * max(0.0, x[j])
    if tmp1 > tmp2:
        return 0.0, max(0.0, A[i][j].upper)
    elif tmp1 == tmp2:
        return 0.5 * min(0.0, A[i][j].lower), 0.5 * max(0.0, A[i][j].lower)
    else:
        return min(0.0, A[i][j].lower), 0.0


def dx(x, pos=True):
    if x > 0.0:
        return 1.0 if pos else 0.0
    elif x == 0.0:
        return 0.5 if pos else -0.5
    else:
        return 0.0 if pos else -1.0
