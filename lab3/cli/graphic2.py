import numpy as np
import matplotlib.pyplot as plt

# 1. Исходные узлы (из C++)
xs = np.array([0.0, 0.5, 1.0, 1.5, 2.0])
ys = np.array([0.0, 0.97943, 1.8415, 2.4975, 2.9093])
X_star = 0.8

# 2. Функция сплайна (воспроизводим логику для точного графика)
# Вычисляем h
h = np.diff(xs)
n = len(h)

# Прямой ход метода прогонки
alpha = np.zeros(n)
beta = np.zeros(n)
c = np.zeros(n + 1)

for i in range(1, n):
    A = h[i - 1]
    B = 2.0 * (h[i - 1] + h[i])
    C = h[i]
    F = 3.0 * ((ys[i + 1] - ys[i]) / h[i] - (ys[i] - ys[i - 1]) / h[i - 1])

    z = B - A * alpha[i - 1]
    alpha[i] = C / z
    beta[i] = (F - A * beta[i - 1]) / z

# Обратный ход
for i in range(n - 1, 0, -1):
    c[i] = beta[i] - alpha[i] * c[i + 1]

# Собираем коэффициенты участков
splines = []
for i in range(n):
    a_i = ys[i]
    b_i = (ys[i + 1] - ys[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0
    c_i = c[i]
    d_i = (c[i + 1] - c[i]) / (3.0 * h[i])
    splines.append((xs[i], xs[i + 1], a_i, b_i, c_i, d_i))

# 3. Построение графика
plt.figure(figsize=(11, 7))

# Цвета для разных участков сплайна, чтобы подчеркнуть кусочность
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']

# Отрисовка каждого сегмента сплайна отдельно
for i, (x_start, x_end, a, b, c_coeff, d) in enumerate(splines):
    x_local = np.linspace(x_start, x_end, 100)
    delta = x_local - x_start
    y_local = a + b * delta + c_coeff * (delta ** 2) + d * (delta ** 3)

    plt.plot(x_local, y_local, label=f'Участок {i+1} $[{x_start}; {x_end}]$',
             color=colors[i % len(colors)], linewidth=2.5)

# 4. Вычисление и нанесение точки X*
# Находим нужный интервал для X*
target_val = 0.0
for x_start, x_end, a, b, c_coeff, d in splines:
    if x_start <= X_star <= x_end:
        delta = X_star - x_start
        target_val = a + b * delta + c_coeff * (delta ** 2) + d * (delta ** 3)
        break

# Рисуем точку X* и проекции
plt.scatter([X_star], [target_val], color='purple', s=160, zorder=6, marker='*',
            label=f'Точка $(X^*, S(X^*)) = ({X_star}, {target_val:.4f})$', edgecolors='black')
plt.vlines(x=X_star, ymin=0, ymax=target_val, color='purple', linestyle=':', alpha=0.6)
plt.hlines(y=target_val, xmin=0, xmax=X_star, color='purple', linestyle=':', alpha=0.6)

# 5. Исходные узлы интерполяции
plt.scatter(xs, ys, color='black', s=80, zorder=5, label='Исходные узлы $(x_i, y_i)$')

# Вертикальные пунктирные линии на границах узлов для наглядности склейки
for x_node in xs:
    plt.axvline(x=x_node, color='gray', linestyle='--', alpha=0.3)

# Оформление
plt.title('Интерполяция кубическим сплайном (естественный сплайн)', fontsize=14, fontweight='bold')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='lower right')
plt.xlim(-0.1, 2.1)
plt.ylim(-0.2, 3.2)

plt.tight_layout()
plt.show()