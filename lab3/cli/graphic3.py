import numpy as np
import matplotlib.pyplot as plt

# 1. Исходные данные из лабораторной работы
xs = np.array([-1.0, 0.0, 1.0, 2.0, 3.0, 4.0])
ys = np.array([-1.8415, 0.0, 1.8415, 2.9093, 3.1411, 3.2432])

# 2. МНК для степеней 1, 2, 3 (получение полиномов)
# np.polyfit возвращает коэффициенты от старшей степени к младшей: [a_m, ..., a_1, a_0]
poly_deg1 = np.poly1d(np.polyfit(xs, ys, 1))
poly_deg2 = np.poly1d(np.polyfit(xs, ys, 2))
poly_deg3 = np.poly1d(np.polyfit(xs, ys, 3))

# Расчет MSE для каждого многочлена
mse1 = np.mean((poly_deg1(xs) - ys) ** 2)
mse2 = np.mean((poly_deg2(xs) - ys) ** 2)
mse3 = np.mean((poly_deg3(xs) - ys) ** 2)

# 3. Точки для отрисовки гладких кривых
x_smooth = np.linspace(min(xs) - 0.5, max(xs) + 0.5, 300)

# 4. Построение графика
plt.figure(figsize=(11, 7))

# Отрисовка исходных узлов
plt.scatter(xs, ys, color='black', s=90, zorder=5,
            label='Исходные данные', edgecolors='white', linewidth=1.5)

# Отрисовка многочленов
plt.plot(x_smooth, poly_deg1(x_smooth), label=f'Многочлен 1-й ст. (MSE = {mse1:.4f})',
         color='#1f77b4', linestyle='--', linewidth=2)
plt.plot(x_smooth, poly_deg2(x_smooth), label=f'Многочлен 2-й ст. (MSE = {mse2:.4f})',
         color='#2ca02c', linestyle='-.', linewidth=2)
plt.plot(x_smooth, poly_deg3(x_smooth), label=f'Многочлен 3-й ст. (MSE = {mse3:.4f})',
         color='#d62728', linestyle='-', linewidth=2.5)

# Дополнительно: пунктирные проекции от точек до лучшего многочлена (показывают ошибки/невязки)
for x_val, y_val in zip(xs, ys):
    y_pred = poly_deg3(x_val)
    plt.vlines(x=x_val, ymin=min(y_val, y_pred), ymax=max(y_val, y_pred),
               color='red', linestyle=':', alpha=0.6)

# Оформление
plt.title('Приближение табличной функции методом наименьших квадратов (МНК)', fontsize=14, fontweight='bold')
plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=11, loc='upper left')
plt.xlim(min(xs) - 0.6, max(xs) + 0.6)

plt.tight_layout()
plt.show()