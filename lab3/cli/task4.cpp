#include <array>
#include <cassert>
#include <iomanip>
#include <iostream>
#include <span>
#include <vector>


struct DifferentiationResult {
    double first_derivative;
    double second_derivative;
    std::vector<double> used_nodes;
};

auto Differentiation(
    const std::span<double> &x_values,
    const std::span<double> &y_values,
    double x
) -> DifferentiationResult {
    assert(x_values.size() == y_values.size());

    const auto n = x_values.size();

    assert(n >= 3);

    auto closest_index = 0ll;
    auto min_distance = std::abs(x_values[0] - x);

    for (std::uint64_t i = 1; i < n; ++i) {
        if (std::abs(x_values[i] - x) < min_distance) {
            min_distance = std::abs(x_values[i] - x);

            closest_index = static_cast<std::int64_t>(i);
        }
    }

    auto start_index = closest_index - 1;

    if (start_index < 0) {
        start_index = 0;
    } else if (start_index + 2 >= n) {
        start_index = static_cast<std::int64_t>(n) - 3;
    }

    const auto x_0 = x_values[start_index],     y_0 = y_values[start_index];
    const auto x_1 = x_values[start_index + 1], y_1 = y_values[start_index + 1];
    const auto x_2 = x_values[start_index + 2], y_2 = y_values[start_index + 2];

    const auto f01 = (y_1 - y_0) / (x_1 - x_0);
    const auto f12 = (y_2 - y_1) / (x_2 - x_1);
    const auto f012 = (f12 - f01) / (x_2 - x_0);

    const auto first_derivative = f01 + f012 * (2 * x - x_0 - x_1);
    const auto second_derivative = 2 * f012;

    return {
        .first_derivative = first_derivative,
        .second_derivative = second_derivative,
        .used_nodes = {x_0, x_1, x_2}
    };
}

int main(int argc, char **argv) {
    // 3.4. Численное дифференцирование
    //
    // Вычислить первую и вторую производную от таблично заданной функции y_i = f(x_i), где i = 0,...,4
    // в точке x = X^{*}.
    //
    // В результате работы программы необходимо вывести:
    // 1. Исходную таблицу значений функции.
    // 2. Точку, в которой выполняется вычисление.
    // 3. Узлы таблицы, использованные при вычислении производных.
    // 4. Значение первой производной.
    // 5. Значение второй производной.

    std::array<double, 5> xs = {
        -1.0, -0.4, 0.2, 0.6, 1.0
    };
    std::array<double, 5> ys = {
        -1.4142, -0.55838, 0.27870, 0.84008, 1.4142
    };
    const auto X = 0.2;

    const auto orders = {1, 2};

    // 1. Исходная таблица значений функции.
    std::cout << "1. Function Table\n\n";
    std::cout << "i" << std::setw(5) << "x" << std::setw(10) << "f(x)" << "\n";
    std::cout << std::string(16, '-') << "\n";
    for (std::uint64_t i = 0; i < xs.size(); ++i) {
        std::cout << i + 1 << std::setw(5) << xs[i] << std::setw(10) << ys[i] << "\n";
    }

    // 2. Точка, в которой выполняется вычисление.
    std::cout << "\n2. X = " << X << "\n";

    auto [first_derivative, second_derivative, used_nodes] = Differentiation(xs, ys, X);

    // 3. Узлы таблицы, использованные при вычислении производных.
    std::cout << "\n3. Used nodes: ";
    for (const auto &used_node : used_nodes) {
        std::cout << used_node << " ";
    }
    std::cout << "\n";

    // 4. Значение первой производной.
    std::cout << "\n4. First derivative: " << first_derivative << "\n";

    // 5. Значение второй производной.
    std::cout << "\n5. Second derivative: " << second_derivative << "\n";

    return 0;
}
