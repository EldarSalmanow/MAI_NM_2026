import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 3**x - 5 * x**2 + 1

x = np.linspace(-2, 4.5, 1000)

plt.figure(figsize=(8, 5))
plt.plot(x, f(x), label=r'$f(x) = 3^x - 5x^2 + 1$', color='blue', lw=2)
plt.axhline(0, color='black', linestyle='--', linewidth=1)  # Ось OX
plt.grid(True, linestyle=':', alpha=0.7)

plt.title('Локализация корней: пересечения с осью OX')
plt.xlabel('x')
plt.ylabel('y')
plt.ylim(-15, 15)
plt.legend()
plt.show()
