#include <array>
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <span>
#include <sstream>
#include <vector>


auto LagrangeEvaluate(
    const std::span<double> &x_values,
    const std::span<double> &y_values,
    double x
) -> double {
    assert(x_values.size() == y_values.size());

    const auto n = x_values.size();

    auto interpolation = 0.0;

    for (std::size_t i = 0; i < n; ++i) {
        auto l_i = 1.0;

        for (std::size_t j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }

            l_i *= (x - x_values[j]) / (x_values[i] - x_values[j]);
        }

        interpolation += y_values[i] * l_i;
    }

    return interpolation;
}

auto LagrangePolynomialString(
    const std::span<double> &x_values,
    const std::span<double> &y_values
) -> std::string {
    assert(x_values.size() == y_values.size());

    const auto n = x_values.size();

    std::stringstream string_stream;

    string_stream << "L(x) = ";

    for (std::size_t i = 0; i < n; ++i) {
        auto denominator = 1.0;

        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) {
                denominator *= (x_values[i] - x_values[j]);
            }
        }

        const auto coefficient = y_values[i] / denominator;

        if (i > 0 && coefficient >= 0) {
            string_stream << " + ";
        }

        string_stream << coefficient;

        for (std::size_t j = 0; j < n; ++j) {
            if (i != j) {
                string_stream << "*(x - " << x_values[j] << ")";
            }
        }
    }

    return string_stream.str();
}

auto DividedDifferences(
    const std::span<double> &x_values,
    const std::span<double> &y_values
) -> std::vector<double> {
    assert(x_values.size() == y_values.size());

    const auto n = x_values.size();

    std::vector<double> differences(y_values.begin(), y_values.end());
    for (std::size_t k = 1; k < n; ++k) {
        for (std::size_t i = n - 1; i >= k; --i) {
            differences[i] = (differences[i] - differences[i - 1]) / (x_values[i] - x_values[i - k]);
        }
    }

    return differences;
}

auto NewtonEvaluate(
    const std::span<double> &x_values,
    const std::span<double> &y_values,
    const std::span<double> &differences,
    double x
) -> double {
    const auto n = x_values.size();

    auto interpolation = differences.back();

    for (std::size_t i = n - 1; i > 0; --i) {
        const auto index = i - 1;

        interpolation = interpolation * (x - x_values[index]) + differences[index];
    }

    return interpolation;
}

auto NewtonPolynomialString(
    const std::span<double> &x_values,
    const std::span<double> &differences
) -> std::string {
    assert(x_values.size() == differences.size());

    const auto n = x_values.size();

    std::stringstream string_stream;

    string_stream << "N(x) = ";

    string_stream << differences[0];
    for (std::size_t i = 1; i < n; ++i) {
        const auto difference = differences[i];

        if (difference >= 0) {
            string_stream << " + ";
        }

        string_stream << difference;

        for (std::size_t j = 0; j < i; ++j) {
            string_stream << "*(x - " << x_values[j] << ")";
        }
    }

    return string_stream.str();
}

int main(int argc, char **argv) {
    // 3.1. Интерполяционные многочлены Лагранжа и Ньютона
    //
    // Используя таблицу значений y_i функции y = f(x), вычисленных в точках x_i, где i = 0,...,3
    // построить интерполяционные многочлены Лагранжа и Ньютона, проходящие через точки {x_i, y_i}.
    // Вычислить значение погрешности интерполяции в точке X^{*}.
    //
    // Необходимо реализовать оба метода.
    // В результате работы программы необходимо вывести:
    // 1. Исходные узлы интерполяции x_i и значения функции y_i.
    // 2. Интерполяционный многочлен Лагранжа.
    // 3. Таблицу разделённых разностей.
    // 4. Интерполяционный многочлен Ньютона.
    // 5. Значение многочлена Лагранжа в заданной точке X^{*}.
    // 6. Значение многочлена Ньютона в заданной точке X^{*}.
    // 7. Точное значение функции f(X^{*}).
    // 8. Абсолютную погрешность интерполяции.
    // 9. Разность между результатами методов Лагранжа и Ньютона.
    //
    // Также необходимо построить общий график, на котором должны быть показаны:
    // - исходная функция;
    // - многочлен Лагранжа;
    // - многочлен Ньютона;
    // - узлы интерполяции;
    // - заданная точка X^{*}.

    const auto function = [] (const double x) -> double {
        return std::sin(x) + x;
    };
    auto lagrangeXs = std::array<double, 4> {
        0.0,
        std::numbers::pi / 6.0,
        std::numbers::pi / 3.0,
        std::numbers::pi / 2.0
    };
    auto newtonXs = std::array<double, 4> {
        0.0,
        std::numbers::pi / 6.0,
        std::numbers::pi / 4.0,
        std::numbers::pi / 2.0
    };
    const auto X = 1.0;

    auto x_nodes = newtonXs;

    std::array<double, 4> y_nodes {};
    for (std::uint64_t i = 0; i < x_nodes.size(); ++i) {
        y_nodes[i] = function(x_nodes[i]);
    }

    std::cout << "1. Nodes\n\n";
    for (std::size_t i = 0; i < x_nodes.size(); ++i) {
        std::cout << "x[" << i << "] = " << std::setw(8) << x_nodes[i] << ", y[" << i << "] = " << std::setw(8) << y_nodes[i] << "\n";
    }

    std::cout << "\n2. Lagrange Polynom: " << LagrangePolynomialString(x_nodes, y_nodes) << "\n";

    auto differences = DividedDifferences(x_nodes, y_nodes);

    std::cout << "\n3. Divided Differences Table\n";
    for (std::uint64_t i = 0; i < differences.size(); ++i) {
        std::cout << "   f(";

        for (std::uint64_t j = 0; j <= i; ++j) {
            std::cout << (j == 0 ? "" : ", ") << "x" << j;
        }

        std::cout << ") = " << differences[i] << "\n";
    }

    std::cout << "\n4. Newton Polynom: " << NewtonPolynomialString(x_nodes, differences) << "\n";

    const auto lagrange = LagrangeEvaluate(x_nodes, y_nodes, X);
    const auto newton = NewtonEvaluate(x_nodes, y_nodes, differences, X);
    const auto y = function(X);

    std::cout << "\n5. L(X) = " << lagrange << "\n";
    std::cout << "\n6. N(X) = " << newton << "\n";
    std::cout << "\n7. f(X) = " << y << "\n";
    std::cout << "\n8. Absolute Error L(X): " << std::abs(y - lagrange) <<
                 "\n   Absolute Error N(X): " << std::abs(y - newton) << "\n";
    std::cout << "\n9. Difference |L(X) - N(X)| = " << std::abs(lagrange - newton) << "\n";

    return 0;
}
