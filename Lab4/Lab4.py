import numpy as np
from numpy.linalg import lstsq

# Функция для вычисления коэффициентов методом наименьших квадратов
def solve_lsm_approximation(x_values, y_values):
    x = np.array(x_values)
    y = np.array(y_values)
    
    # Приведение уравнения к линейной форме: P = a + b * v^2
    X = np.vstack([np.ones_like(x), x]).T  # Матрица [1, x]
    
    # Решение методом наименьших квадратов
    coeffs, _, _, _ = lstsq(X, y, rcond=None)
    return coeffs

# Функция для вычисления стандартного отклонения
def calculate_standard_deviation(x_values, y_values, coefficients):
    x = np.array(x_values)
    y = np.array(y_values)
    
    # Прогнозируемые значения
    predicted = coefficients[0] + coefficients[1] * x
    
    # Вычисление остаточной суммы квадратов
    residual_sum = np.sum((y - predicted)**2)
    n = len(y)
    degree = len(coefficients) - 1
    
    # Стандартное отклонение
    return np.sqrt(residual_sum / (n - degree - 1))

if __name__ == "__main__":
    experimental_x = [2.40**2, 3.50**2, 5.00**2, 6.89**2, 10.00**2]  # x = v^2
    experimental_y = [0.0141, 0.0281, 0.0562, 0.1125, 0.2250]  # P значения

    # Решение задачи методом наименьших квадратов
    coefficients = solve_lsm_approximation(experimental_x, experimental_y)

    print("The equation: P = a + b * v^2")
    print(f"The coefficients found: a = {coefficients[0]}, b = {coefficients[1]}")

    # Вычисление стандартного отклонения
    std_deviation = calculate_standard_deviation(experimental_x, experimental_y, coefficients)
    print(f"\nThe standard deviation of the equation: {std_deviation}")
