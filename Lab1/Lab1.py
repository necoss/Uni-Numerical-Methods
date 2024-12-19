import numpy as np

def gaussElimination(A, b):
    n = len(b)
    A = A.astype(float)
    b = b.astype(float)

    for k in range(n):
        # Выбор главного элемента (строка с макс модулем)
        max_row = np.argmax(np.abs(A[k:, k])) + k # Находим наибольшее абсолютное значение
        A[[k, max_row]] = A[[max_row, k]]
        b[[k, max_row]] = b[[max_row, k]] # Меняем для устойчивости вычислений

        # Обнуление элементов ниже главного элемента
        for i in range(k + 1, n):
            factor = A[i, k] / A[k, k] # для строк где i > k мы вычисляем коэфф factor
            A[i, k:] -= factor * A[k, k:]
            b[i] -= factor * b[k]

    # Обратный ход (формула)
    x = np.zeros(n)
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - np.dot(A[i, i + 1:], x[i + 1:])) / A[i, i]

    return x

def ldlFactorization(A): # Разложить симметричную матрицу A в произведение LDL^T
    n = A.shape[0]
    L = np.zeros_like(A, dtype=float)
    D = np.zeros(n, dtype=float)

    for i in range(n):
        for j in range(i):
            L[i, j] = A[i, j] - np.sum(L[i, :j] * D[:j] * L[j, :j])
            L[i, j] /= D[j]
        L[i, i] = 1  # Диагональные элементы L = 1
        D[i] = A[i, i] - np.sum(L[i, :i]**2 * D[:i])

    return L, D

def solveLdl(L, D, b):
    # нижнетреугольная матрица
    y = np.zeros_like(b, dtype=float)
    for i in range(len(b)):
        y[i] = b[i] - np.dot(L[i, :i], y[:i])

    # Solve Dz = y
    z = y / D

    # Обратный ход L^T x = z
    x = np.zeros_like(b, dtype=float)
    for i in range(len(b) - 1, -1, -1):
        x[i] = z[i] - np.dot(L[i + 1:, i], x[i + 1:])

    return x


A = np.array([
    [8.64, -6.39, 4.21],
    [1.71, 4.25, 7.92],
    [5.42, 1.84, -3.41]
])
b = np.array([10.21, 3.41, 12.29])

x_gauss = gaussElimination(A.copy(), b.copy())
print("Solution using Gaussian Elimination:", x_gauss)

L, D = ldlFactorization(A)
x_ldl = solveLdl(L, D, b)
print("Solution using LDL^T Factorization:", x_ldl)