import numpy as np

def jacobian(x):
    """
    Вычисляет матрицу Якоби для системы нелинейных уравнений.
    """
    x1, x2 = x
    J = np.array([
        [6 * x1**2, -1],
        [2 * x1, 2 * x2]
    ])
    return J

def equations(x):
    """
    Вычисляет значения системы уравнений.
    """
    x1, x2 = x
    return np.array([
        2 * x1**3 - x2 - 5,
        x1**2 + x2**2 - 9
    ])

def newton_method(f, J, x0, tol=1e-6, max_iter=100):
    """
    Метод Ньютона для решения системы нелинейных уравнений.

    f: функция, возвращающая значения системы уравнений.
    J: функция, возвращающая матрицу Якоби.
    x0: начальное приближение.
    tol: допустимая погрешность.
    max_iter: максимальное количество итераций.

    Возвращает решение или сообщение об ошибке.
    """
    x = x0
    for i in range(max_iter):
        Fx = f(x)
        Jx = J(x)

        # Решаем линейную систему Jx * delta_x = -Fx
        delta_x = np.linalg.solve(Jx, -Fx)

        # Обновляем приближение
        x = x + delta_x

        # Проверяем критерий сходимости
        if np.linalg.norm(delta_x, ord=np.inf) < tol:
            return x, i + 1

    raise ValueError("Newton's method did not converge for {} iterations".format(max_iter))

# Начальное приближение
x0 = np.array([1.0, 1.0])

# Решение методом Ньютона
try:
    solution, iterations = newton_method(equations, jacobian, x0)
    print("The solution was found for {} iterations: {}".format(iterations, solution))
except ValueError as e:
    print(e)
