import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange

# 1. Исходные данные (те же, что и в C++ коде)
x_nodes = np.array([0.0, np.pi/6.0, np.pi/4.0, np.pi/2.0])
y_nodes = np.sin(x_nodes) + x_nodes
X_star = 1.0

def exact_function(x):
    return np.sin(x) + x

# 2. Создаем массив точек для отрисовки гладких кривых
# Берем интервал чуть шире, чем крайние узлы, чтобы было видно поведение
x_smooth = np.linspace(-0.2, 2.0, 300)
y_exact = exact_function(x_smooth)

# 3. Вычисление многочленов
# Используем встроенную функцию Лагранжа из SciPy (она математически идентична твоей из C++)
poly = lagrange(x_nodes, y_nodes)
y_lagrange = poly(x_smooth)

# Многочлен Ньютона математически равен Лагранжу, поэтому кривая будет той же
# newton(x_nodes, y_nodes)
y_newton = newton(x_nodes, y_nodes)

# 4. Настройка графика
plt.figure(figsize=(10, 6))

# Рисуем исходную функцию
plt.plot(x_smooth, y_exact, label='Исходная функция $f(x) = \sin(x) + x$',
         color='black', linewidth=2)

# Рисуем многочлен Лагранжа (толстая полупрозрачная синяя линия)
plt.plot(x_smooth, y_lagrange, label='Многочлен Лагранжа',
         color='blue', linewidth=4, alpha=0.5)

# Рисуем многочлен Ньютона (красный пунктир поверх Лагранжа, чтобы показать совпадение)
plt.plot(x_smooth, y_newton, label='Многочлен Ньютона',
         color='red', linestyle='--', linewidth=2)

# Отмечаем узлы интерполяции
plt.scatter(x_nodes, y_nodes, color='green', s=100, zorder=5,
            label='Узлы интерполяции', edgecolors='black')

# Отмечаем заданную точку X* (и ее значение на интерполянте)
Y_star_calc = poly(X_star)
plt.scatter([X_star], [Y_star_calc], color='magenta', s=150, zorder=5, marker='*',
            label=f'Точка $X^* = {X_star}$', edgecolors='black')

# Вертикальная линия для наглядности точки X*
plt.vlines(x=X_star, ymin=0, ymax=Y_star_calc, color='magenta', linestyle=':', alpha=0.7)

# Оформление
plt.title('Интерполяция многочленами Лагранжа и Ньютона', fontsize=14, fontweight='bold')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10, loc='upper left')
plt.xlim(-0.2, 2.0)
plt.ylim(0, 3.0)

# Показываем график
plt.tight_layout()
plt.show()