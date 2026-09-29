import numpy as np
import matplotlib.pyplot as plt

# Параметр a = 3
t = np.linspace(0, 2 * np.pi, 500)

# 1-е уравнение: x1 = 3 + cos(x2)
x1_eq1 = 3 + np.cos(t)
x2_eq1 = t

# 2-е уравнение: x2 = 3 + sin(x1)
x1_eq2 = t
x2_eq2 = 3 + np.sin(t)

plt.figure(figsize=(6, 6))
plt.plot(x1_eq1, x2_eq1, label=r'$x_1 - \cos(x_2) = 3$', color='purple', lw=2)
plt.plot(x1_eq2, x2_eq2, label=r'$x_2 - \sin(x_1) = 3$', color='orange', lw=2)

plt.grid(True, linestyle=':', alpha=0.7)
plt.title('Точка пересечения — решение системы')
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')
plt.xlim(1.5, 4.5)
plt.ylim(1.5, 4.5)
plt.legend()
plt.show()
