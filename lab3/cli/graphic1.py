import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import lagrange#, newton

x_nodes = np.array([0.0, np.pi/6.0, np.pi/3.0, np.pi/2.0])
# x_nodes_2 = np.array([0.0, np.pi/6.0, np.pi/4.0, np.pi/2.0])
y_nodes = np.sin(x_nodes) + x_nodes
X_star = 1.0

def exact_function(x):
    return np.sin(x) + x

x_smooth = np.linspace(-0.2, 2.0, 300)
y_exact = exact_function(x_smooth)

poly = lagrange(x_nodes, y_nodes)
y_lagrange = poly(x_smooth)

# y_newton = newton(x_nodes, y_nodes)

plt.figure(figsize=(10, 6))

plt.plot(x_smooth, y_exact, label='Исходная функция $f(x) = \sin(x) + x$',
         color='black', linewidth=2)

plt.plot(x_smooth, y_lagrange, label='Многочлен Лагранжа',
         color='blue', linewidth=4, alpha=0.5)

# plt.plot(x_smooth, y_newton, label='Многочлен Ньютона',
#          color='red', linestyle='--', linewidth=2)

plt.scatter(x_nodes, y_nodes, color='green', s=100, zorder=5,
            label='Узлы интерполяции', edgecolors='black')

Y_star_calc = poly(X_star)
plt.scatter([X_star], [Y_star_calc], color='magenta', s=150, zorder=5, marker='*',
            label=f'Точка $X^* = {X_star}$', edgecolors='black')

plt.vlines(x=X_star, ymin=0, ymax=Y_star_calc, color='magenta', linestyle=':', alpha=0.7)

plt.title('Интерполяция многочленами Лагранжа и Ньютона', fontsize=14, fontweight='bold')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(fontsize=10, loc='upper left')
plt.xlim(-0.2, 2.0)
plt.ylim(0, 3.0)

plt.tight_layout()
plt.show()