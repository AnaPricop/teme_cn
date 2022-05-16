import numpy as np

eps = 1e-8


def read_files(file_a, file_b):
    diag = []
    with open(file_a) as f:
        size = int(f.readline())
        next(f)
        a = [[] for _ in range(size)]

        for line in f:
            line_sp = line.split(', ')
            try:
                value = float(line_sp[0])
            except:
                return diag, a, size
            i = int(line_sp[1])
            j = int(line_sp[2])
            if validare_i_j((i, j)):
                if if_exists_i_j((i, j), a):  # daca gasim aceeasi indexi, adaugam la val initiala
                    i2 = i
                    j2 = j
                    line = a[i2]
                    for element in line:
                        if element[1] == j2:
                            element[0] += value
                else:
                    entry = (value, j)
                    a[i].append(entry)
            if i == j:  # daca gasim element de pe diagonala
                diag.append(value)
    with open(file_b) as f:
        b = []
        for line in f:
            try:
                value = float(line)
                b.append(value)
            except:
                return b
    return b, a, diag, size


def find_ij(line, p):
    for element in line:
        if element[1] == p:
            return element
    return None


def if_exists_i_j(p, a):
    i = p[0]
    j = p[1]
    if find_ij(a[i], j) != None:
        return True
    return False


def validare_i_j(p):
    if p[0] < p[1]:
        print('am gasit j>i')
        return False
    return True


def verificare_diagonala_a(a, len):
    for i in range(len):
        if find_ij(a[i], i) == None:
            return False
    return True


mat_b1, mat_a1, d1, len1 = read_files('a_1.txt', 'b_1.txt')

mat_b2, mat_a2, d2, len2 = read_files('a_2.txt', 'b_2.txt')

mat_b3, mat_a3, d3, len3 = read_files('a_3.txt', 'b_3.txt')

mat_b4, mat_a4, d4, len4 = read_files('a_4.txt', 'b_4.txt')

mat_b5, mat_a5, d5, len5 = read_files('a_5.txt', 'b_5.txt')

if not verificare_diagonala_a(mat_a1, len1):
    print('Am gasit element nul pe diagonala matricei a, deci nu putem aplica Jacobi.')
if not verificare_diagonala_a(mat_a2, len2):
    print('Am gasit element nul pe diagonala matricei a, deci nu putem aplica Jacobi.')
if not verificare_diagonala_a(mat_a3, len3):
    print('Am gasit element nul pe diagonala matricei a3, deci nu putem aplica Jacobi.')
if not verificare_diagonala_a(mat_a4, len4):
    print('Am gasit element nul pe diagonala matricei a4, deci nu putem aplica Jacobi.')
if not verificare_diagonala_a(mat_a5, len5):
    print('Am gasit element nul pe diagonala matricei a5, deci nu putem aplica Jacobi.')


def dot_product(a, x):
    result = [0 for _ in range(len(x))]
    for i in range(len(a)):
        line = a[i]
        for aux in line:
            val = aux[0]
            j = aux[1]
            result[i] += val * x[j]
            if i != j:  # not on diag
                result[j] += val * x[i]
    return result


def verify_solution(a, b, x):
    dot = dot_product(a, x)
    return abs(max(np.subtract(dot, b), key=abs))


def jacobi_method(a, b, d, len):
    x_c = [0 for _ in range(len)] # pseudocod
    x_p = [0 for _ in range(len)]
    k = 0
    delta_x = eps
    k_max = 10000
    while delta_x >= eps and k <= k_max and delta_x <= 1e+8:
        x_p = x_c
        x_c = jacobi_st(x_p, a, b, d, len)
        delta_x = np.linalg.norm(np.subtract(x_c, x_p))
        k = k + 1
    if delta_x < eps:
        print('Iterations:', k)
        return x_c
    if not verificare_diagonala_a(a, len):
        return None
    return None


def jacobi_st(x_prev, a, b, d, size):
    x = [b[i] for i in range(size)]
    for i in range(size):
        for element in a[i]:
            value = element[0]
            j = element[1]
            if i != j:
                x[i] -= value * x_prev[j]
                x[j] -= value * x_prev[i]
    for i in range(size):
        x[i] /= d[i]

    return x


xj1 = jacobi_method(mat_a1, mat_b1, d1, len1)
if xj1 is None:
    print('divergence')
else:
    norma1 = verify_solution(mat_a1, mat_b1, xj1)
    print('Norma 1:', norma1)

xj2 = jacobi_method(mat_a2, mat_b2, d2, len2)
if xj2 is None:
    print('divergence')
else:
    norma2 = verify_solution(mat_a2, mat_b2, xj2)
    print('Norma 2:', norma2)

xj3 = jacobi_method(mat_a3, mat_b3, d3, len3)
if xj3 is None:
    print('divergence')
else:
    norma3 = verify_solution(mat_a3, mat_b3, xj3)
    print('Norma 3:', norma3)

xj4 = jacobi_method(mat_a4, mat_b4, d4, len4)
if xj4 is None:
    print('divergence')
else:
    norma4 = verify_solution(mat_a4, mat_b4, xj4)
    print('Norma 4:', norma4)
#
# xj5 = jacobi_method(mat_a5, mat_b5, d5, len5)
# if xj5 is None:
#     print('divergence')
# else:
#     norma5 = verify_solution(mat_a5, mat_b5, xj5)
#     print('Norma 5:', norma5)
