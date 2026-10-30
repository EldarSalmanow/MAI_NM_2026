#include <array>
#include <iomanip>
#include <iostream>
#include <span>
#include <tuple>
#include <vector>


struct SplineSegment {
    double x_1, x_2;
    double a, b, c, d;
};

auto CubicSpline(
    const std::span<double> x_values,
    const std::span<double> y_values
) -> std::vector<SplineSegment> {
    const auto n = x_values.size() - 1;

    std::vector<double> h(n);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = x_values[i + 1] - x_values[i];
    }

    std::vector<double> c(n + 1, 0.0);
    std::vector<double> p(n, 0.0);
    std::vector<double> q(n, 0.0);

    for (std::size_t i = 1; i < n; ++i) {
        const double A = h[i - 1];
        const double B = 2.0 * (h[i - 1] + h[i]);
        const double C = h[i];
        const double F = 3.0 * ((y_values[i + 1] - y_values[i]) / h[i] - (y_values[i] - y_values[i - 1]) / h[i - 1]);

        const double z = B - A * p[i - 1];
        p[i] = C / z;
        q[i] = (F - A * q[i - 1]) / z;
    }

    for (std::size_t i = n - 1; i > 0; --i) {
        c[i] = q[i] - p[i] * c[i + 1];
    }

    std::vector<SplineSegment> spline(n);
    for (std::size_t i = 0; i < n; ++i) {
        spline[i] = SplineSegment {
            .x_1 = x_values[i],
            .x_2 = x_values[i + 1],
            .a = y_values[i],
            .b = (y_values[i + 1] - y_values[i]) / h[i] - h[i] * (2.0 * c[i] + c[i + 1]) / 3.0,
            .c = c[i],
            .d = (c[i + 1] - c[i]) / (3.0 * h[i])
        };
    }

    return spline;
}

auto EvalCubicSpline(
    const std::span<SplineSegment> spline,
    const double x
) -> std::tuple<double, std::uint64_t> {
    std::size_t index = 0;

    if (x <= spline.front().x_1) {
        index = 0;
    } else if (x >= spline.back().x_2) {
        index = spline.size() - 1;
    } else {
        for (std::size_t i = 0; i < spline.size(); ++i) {
            if (spline[i].x_1 <= x && x <= spline[i].x_2) {
                index = i;

                break;
            }
        }
    }

    const auto [x_1, x_2, a, b, c, d] = spline[index];

    const auto delta = x - x_1;

    const auto result = a + b * delta + c * delta * delta + d * delta * delta * delta;

    return {result, index};
}

auto SplineSecondDerivative(
    const SplineSegment &segment,
    double x
) -> double {
    // S''(x) = 2c + 6d(x - x_1)

    return 2.0 * segment.c + 6.0 * segment.d * (x - segment.x_1);
}

int main(int argc, char **argv) {
    // 3.2. Кубический сплайн
    //
    // Построить кубический сплайн для функции, заданной в узлах интерполяции, предполагая, что сплайн
    // имеет нулевую кривизну при x = x_0 и x = x_4. Вычислить значение функции в точке x = X^{*}.
    //
    // В результате работы программы необходимо вывести:
    // 1. Исходные узлы и значения функции.
    // 2. Коэффициенты a_i, b_i, c_i, d_i для каждого участка сплайна.
    // 3. Интервал, на котором действует каждый полученный многочлен.
    // 4. Значение сплайна в заданной точке X^{*}.
    // 5. Номер интервала, содержащего точку X^{*}.
    // 6. Проверку прохождения сплайна через исходные узлы.
    // 7. Значения второй производной на концах интервала для проверки граничных условий.
    //
    // Также необходимо построить график, на котором должны быть показаны:
    // - исходные узлы;
    // - все участки кубического сплайна;
    // - заданная точка (X^{*}, S(X^{*})).

    std::array<double, 5> xs = {
        0.0, 0.5, 1.0, 1.5, 2.0
    };
    std::array<double, 5> ys = {
        0.0, 0.97943, 1.8415, 2.4975, 2.9093
    };
    const double X = 0.8;

    // 1. Исходные узлы и значения функции.
    std::cout << "1. Nodes\n\n";
    for (std::size_t i = 0; i < xs.size(); ++i) {
        std::cout << "x[" << i << "] = " << std::setw(8) << xs[i] << ", y[" << i << "] = " << std::setw(8) << ys[i] << "\n";
    }

    auto spline = CubicSpline(xs, ys);

    // 2. Коэффициенты a_i, b_i, c_i, d_i для каждого участка сплайна.
    // 3. Интервал, на котором действует каждый полученный многочлен.
    std::cout << "\n2. & 3. Coefficients & Interval Table\n";
    std::cout << "i" << std::setw(15) << "[x_1; x_2]"
              << std::setw(9) << "a" << std::setw(11) << "b"
              << std::setw(10) << "c" << std::setw(12) << "d" << "\n";
    std::cout << std::string(63, '-') << "\n";
    for (std::uint64_t i = 0; i < spline.size(); ++i) {
        const auto [x_1, x_2, a, b, c, d] = spline[i];

        std::cout << i + 1 << std::setw(7) << "[" << x_1 << "; " << x_2 << "]"
                  << std::setw(12) << a << std::setw(12) << b
                  << std::setw(12) << c << std::setw(12) << d << "\n";
    }

    const auto [value, interval] = EvalCubicSpline(spline, X);

    // 4. Значение сплайна в заданной точке X^{*}.
    std::cout << "\n4. S(X^{*}) = " << value << "\n";

    // 5. Номер интервала, содержащего точку X^{*}.
    std::cout << "\n5. i = " << interval + 1 << "\n";

    // 6. Проверку прохождения сплайна через исходные узлы.
    std::cout << "\n6. Check spline\n";
    for (std::size_t i = 0; i < xs.size(); ++i) {
        const auto [val, _] = EvalCubicSpline(spline, xs[i]);
        std::cout << "S(" << xs[i] << ") = " << val << " (expected " << ys[i] << "), "
                  << "error = " << std::abs(val - ys[i]) << "\n";
    }

    // 7. Значения второй производной на концах интервала для проверки граничных условий.
    std::cout << "\n7. Check borders\n";
    if (!spline.empty()) {
        std::cout << "S''(x_0) = " << SplineSecondDerivative(spline.front(), spline.front().x_1) << "\n";
        std::cout << "S''(x_n) = " << SplineSecondDerivative(spline.back(), spline.back().x_2) << "\n";
    }

    return 0;
}
