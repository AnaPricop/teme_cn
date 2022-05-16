import sys

import numpy as np


def read_file():
    A = []
    with open('matrix_1.txt', 'r') as rd:
        ln = rd.readline()
        while ln:
            A.append([float(element) for element in ln.split(' ')])
            ln = rd.readline()
    A = np.array(A)
    A_init = A.copy()
    return A, A_init


def index(A):
    n = len(A)
    max = 0
    p = -1
    q = -1
    for j in range(n):
        for i in range(j):
            if max < abs(A[i][j]):
                max = abs(A[i][j])
                p = i
                q = j
    return p, q


def check_value(value, eps):
    if abs(value) < abs(eps):  # verificam daca valoarea e aproape de 0
        return 0

    return value


def jacobi_algorithm(n, A, k_max, eps):
    # pseudocod tema
    k = 0
    p, q = index(A)
    U = np.identity(n)

    alpha = (A[p][p] - A[q][q]) / (2 * A[p][q])
    t = -alpha + (1 if alpha >= 0 else -1) * ((alpha ** 2 + 1) ** (1 / 2))
    c = 1 / ((1 + t ** 2) ** (1 / 2))
    s = t / ((1 + t ** 2) ** (1 / 2))

    while not este_diagonala(A, eps) and k <= k_max:

        for j in range(n):
            if j != p and j != q:
                A[p][j] = check_value(c * A[p][j] + s * A[q][j], eps)

        for j in range(n):
            if j != p and j != q:
                A[q][j] = A[j][q] = check_value(-s * A[j][p] + c * A[q][j], eps)

        for j in range(n):
            if j != p and j != q:
                A[j][p] = check_value(A[p][j], eps)

        A[p][p] = check_value(A[p][p] + t * A[p][q], eps)
        A[q][q] = check_value(A[q][q] - t * A[p][q], eps)
        A[p][q] = 0
        A[q][p] = 0

        U_orig = np.copy(U)
        for i in range(n):
            U[i][p] = check_value(c * U[i][p] + s * U[i][q], eps)
            U[i][q] = check_value(-s * U_orig[i][p] + c * U[i][q], eps)

        p, q = index(A)

        if p == -1 and q == -1:
            break

        alpha = (A[p][p] - A[q][q]) / (2 * A[p][q])
        t = -alpha + (1 if alpha >= 0 else -1) * ((alpha ** 2 + 1) ** (1 / 2))
        c = 1 / ((1 + t ** 2) ** (1 / 2))
        s = t / ((1 + t ** 2) ** (1 / 2))
        k += 1

    return U


def verifica_norma(jac, lib):
    s = 0
    # elementele calculate cu jacobi
    for diag_jacobi in jac:
        minimum = sys.maxsize
        # elem calculate prin librarie
        for diag_lib in lib:
            if abs(diag_jacobi - diag_lib) < minimum:
                minimum = abs(diag_jacobi - diag_lib)
        s += minimum
    return s


def este_diagonala(matrix, eps):
    # verif daca matricea e diagonala
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if i != j and abs(matrix[i][j]) > eps:
                return False
    return True


eps = 10 ** -8
k_max = 100
A, A_init = read_file()
n = len(A)
# print(A)
U = jacobi_algorithm(n, A, k_max, eps)
print("Primul subpunct:", end='\n')
print("Valori proprii: ", np.diag(A), end='\n')
print("U: ", U, end='\n')
print("Verificarea Ainit*U = Lambda*A: ", check_value(np.linalg.norm(np.matmul(A_init, U) - np.matmul(U, A)), eps))

print("Al doilea subpunct:", end='\n')
eigen_values, eigen_vectors = np.linalg.eigh(A_init)
print("Valorile proprii (eigenvalues): ", eigen_values, end='\n')
print("Vectorii proprii (eigenvectors): ", '\n', eigen_vectors, end='\n')
print("Suma calculata (val proprii Jacobi si val proprii biblioteca): ",
      check_value(verifica_norma(np.diag(A), eigen_values), eps))

print("Al treilea subpunct:", end='\n')
# np linalg svd - pentru Singular Value Decomposition
u, s, v_t = np.linalg.svd(A_init)
s_i = np.array([[0 if i != j else (0 if check_value(s[j], eps) == 0 else 1 / s[j]) for i in range(len(s))] for j in
                range(len(s))])
A_I = np.matmul(np.matmul(np.transpose(v_t), s_i), np.transpose(u))
# Moore-Penrose cu numpy.linalg.pinv
A_J = np.matmul(np.linalg.pinv(np.matmul(np.transpose(A_init), A_init)), np.transpose(A_init))
print("Valorile singulare ale matricei A: ", s)
print("Rangul matricei A: ", np.linalg.matrix_rank(A))
print("Rangul matricei A: ", np.count_nonzero(s[abs(s) >= eps]))
print("Numarul de conditionare al matricei: ", '\n', A_I)
print("Pseudoinversa Moore-Penrose: ", '\n', A_J)
print("Norma : ", check_value(np.linalg.norm(A_I - A_J), eps))
