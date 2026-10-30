#include <array>
#include <cassert>
#include <iostream>
#include <ranges>
#include <span>
#include <sstream>
#include <vector>


auto LeastSquaresMethod(
    std::span<double> x_values,
    std::span<double> y_values,
    std::uint64_t power
) -> std::vector<double> {
    const auto n = power + 1;

    std::vector<std::vector<double>> coefficients(n, std::vector<double>(n + 1, 0.0));
    for (std::uint64_t i = 0; i < n; ++i) {
        for (std::uint64_t j = 0; j < n; ++j) {
            auto temp_a = 0.0;

            for (std::uint64_t k = 0; k < x_values.size(); ++k) {
                temp_a += std::pow(x_values[k], i + j);
            }

            coefficients[i][j] = temp_a;
        }

        auto temp_b = 0.0;

        for (std::uint64_t j = 0; j < x_values.size(); ++j) {
            temp_b += y_values[j] * std::pow(x_values[j], i);
        }

        coefficients[i][n] = temp_b;
    }

    // solve SLE
    for (std::uint64_t i = 0; i < n; ++i) {
        std::uint64_t max_row = i;
        double max_val = std::abs(coefficients[i][i]);

        for (std::uint64_t k = i + 1; k < n; ++k) {
            if (std::abs(coefficients[k][i]) > max_val) {
                max_val = std::abs(coefficients[k][i]);

                max_row = k;
            }
        }

        if (max_row != i) {
            std::swap(coefficients[i], coefficients[max_row]);
        }

        for (std::uint64_t k = i + 1; k < n; ++k) {
            const double factor = coefficients[k][i] / coefficients[i][i];

            for (std::uint64_t j = i; j <= n; ++j) {
                coefficients[k][j] -= factor * coefficients[i][j];
            }
        }
    }

    std::vector<double> result(n, 0.0);
    for (int i = static_cast<int>(n) - 1; i >= 0; --i) {
        result[i] = coefficients[i][n];

        for (std::uint64_t j = i + 1; j < n; ++j) {
            result[i] -= coefficients[i][j] * result[j];
        }

        result[i] /= coefficients[i][i];
    }

    return result;
}

auto MSE(
    const std::span<double> x_values,
    const std::span<double> y_values,
    const std::span<double> coefficients
) -> double {
    assert(x_values.size() == y_values.size());

    const auto n = x_values.size();

    std::vector<double> new_ys(n, 0.0);
    for (std::uint64_t i = 0; i < n; ++i) {
        auto temp_y = 0.0;

        for (std::uint64_t j = 0; j < coefficients.size(); ++j) {
            temp_y += coefficients[j] * std::pow(x_values[i], j);
        }

        new_ys[i] = temp_y;
    }

    auto error = 0.0;

    for (std::uint64_t i = 0; i < new_ys.size(); ++i) {
        error += std::pow(new_ys[i] - y_values[i], 2);
    }

    return error / static_cast<double>(n);
}

auto PolynomialString(
    const std::uint64_t power,
    const std::span<double> coefficients
) -> std::string {
    std::stringstream string_stream;

    string_stream << "P_" << power << "(x) = ";

    for (std::uint64_t i = 0; i < coefficients.size(); ++i) {
        string_stream << (i > 0 ? " + " : "") << coefficients[i] << "x^" << i;
    }

    return string_stream.str();
}

int main(int argc, char **argv) {
    // 3.3. Метод наименьших квадратов
    //
    // Для таблично заданной функции путем решения нормальной системы МНК
    // найти приближающие многочлены a) 1-ой и б) 2-ой степени.
    // Для каждого из приближающих многочленов вычислить сумму квадратов ошибок.
    // Построить графики приближаемой функции и приближающих многочленов.
    //
    // Необходимо построить приближающие многочлены первой, второй и третьей степени.
    // В результате работы программы необходимо вывести:
    // 1. Коэффициенты многочлена первой степени.
    // 2. Коэффициенты многочлена второй степени.
    // 3. Коэффициенты многочлена третьей степени.
    // 4. Полученные многочлены P_1(x), P_2(x) и P_3(x) в явном виде.
    // 5. Значения MSE1, MSE2 и MSE3, рассчитанные по всем исходным точкам.
    // 6. Сравнение полученных ошибок и указание многочлена, который лучше приближает исходные данные.
    //
    // Также необходимо построить общий график, на котором должны быть показаны:
    // - исходные данные в виде отдельных точек;
    // - многочлен первой степени;
    // - многочлен второй степени;
    // - многочлен третьей степени.

    std::array<double, 6> xs {
        -1.0, 0.0, 1.0, 2.0, 3.0, 4.0
    };
    std::array<double, 6> ys {
        -1.8415, 0.0, 1.8415, 2.9093, 3.1411, 3.2432
    };
    const auto powers = {1, 2, 3};

    for (const auto &power : powers) {
        std::cout << power << " power LSM" << "\n";

        auto coefficients = LeastSquaresMethod(xs, ys, power);

        std::cout << "Coefficients: ";
        for (std::uint64_t i = 0; i < coefficients.size(); ++i) {
            std::cout << "a_" << i << " = " << coefficients[i] << "; ";
        }
        std::cout << "\n";

        std::cout << "Polynom: " << PolynomialString(power, coefficients) << "\n";

        std::cout << "MSE: " << MSE(xs, ys, coefficients) << "\n";

        std::cout << "\n";
    }

    return 0;
}