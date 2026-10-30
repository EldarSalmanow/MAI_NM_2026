#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>


template<typename FunctionT>
auto RectangleLeftMethod(
    FunctionT &&function,
    const double a,
    const double b,
    const double h
) -> double {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    auto result = 0.0;

    for (std::uint64_t i = 0; i < n; ++i) {
        const auto x = a + static_cast<double>(i) * h;

        result += function(x);
    }

    return result * h;
}

template<typename FunctionT>
auto RectangleMiddleMethod(
    FunctionT &&function,
    const double a,
    const double b,
    const double h
) -> double {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    auto result = 0.0;

    for (std::uint64_t i = 0; i < n; ++i) {
        const auto x = a + (static_cast<double>(i) + 0.5) * h;

        result += function(x);
    }

    return result * h;
}

template<typename FunctionT>
auto RectangleRightMethod(
    FunctionT &&function,
    const double a,
    const double b,
    const double h
) -> double {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    auto result = 0.0;

    for (std::uint64_t i = 0; i < n; ++i) {
        const auto x = a + static_cast<double>(i + 1) * h;

        result += function(x);
    }

    return result * h;
}

template<typename FunctionT>
auto TrapezoidMethod(
    FunctionT &&function,
    const double a,
    const double b,
    const double h
) -> double {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    auto result = (function(a) + function(b)) / 2.0;

    for (std::uint64_t i = 1; i < n; ++i) {
        const double x = a + static_cast<double>(i) * h;

        result += function(x);
    }

    result *= h;

    return result;
}

template<typename FunctionT>
auto SimpsonMethod(
    FunctionT &&function,
    const double a,
    const double b,
    const double h
) -> double {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    if (n % 2 != 0) {
        std::cout << "[ERROR] Simpson's method requires an even number of intervals!" << std::endl;

        return 0.0;
    }

    auto result = function(a) + function(b);

    for (std::size_t i = 1; i < n; ++i) {
        const double x = a + static_cast<double>(i) * h;

        result += (i % 2 != 0 ? 4.0 : 2.0) * function(x);
    }

    result *= h / 3.0;

    return result;
}

using DefaultFunctionType = double (*) (const double);

template<typename FunctionT>
using MethodFuncType = double (*) (FunctionT &&, const double, const double, const double);

struct MethodData {
    std::string name;
    std::uint32_t p;
    MethodFuncType<DefaultFunctionType> func;
};

auto RungeRombergError(const double I_h1, const double I_h2, const std::uint32_t p) -> double {
    return std::abs(I_h2 - I_h1) / (std::pow(2.0, p) - 1.0);
}

auto RungeRombergRefined(const double I_h1, const double I_h2, const std::uint32_t p) -> double {
    return I_h2 + (I_h2 - I_h1) / (std::pow(2.0, p) - 1.0);
}

int main(int argc, char **argv) {
    // 3.5. Численное интегрирование
    //
    // Вычислить определенный интеграл F = int_{X_0}^{X_1}ydx методами прямоугольников, трапеций, Симпсона с шагами h_1, h_2.
    // Оценить погрешность вычислений, используя Метод Рунге-Ромберга.
    //
    // Необходимо реализовать:
    // * метод левых прямоугольников;
    // * метод правых прямоугольников;
    // * метод средних прямоугольников;
    // * метод трапеций;
    // * метод Симпсона;
    // * уточнение результатов методом Рунге–Ромберга.
    //
    // Результаты необходимо вывести в виде таблицы.
    // Для каждого метода необходимо указать:
    // 1. Значение интеграла с шагом h_1.
    // 2. Значение интеграла с шагом h_2.
    // 3. Порядок точности метода p.
    // 4. Уточнённое значение интеграла по методу Рунге–Ромберга.
    // 5. Оценку погрешности результата, полученного с шагом h_2.
    //
    // Если интеграл можно вычислить аналитически, дополнительно необходимо вывести:
    // - Точное значение интеграла.
    // - Абсолютную погрешность каждого численного метода.
    // - Метод, давший наиболее точный результат.
    //
    // График для этой части не требуется.

    const auto function = [] (const double x) -> double {
        return x / (std::pow(x, 3) + 8.0);
    };

    const auto exact_integral = [] (const double x) -> double {
        return (1.0 / 12.0) * std::log((x * x - 2.0 * x + 4.0) / std::pow(x + 2.0, 2.0)) +
               (std::sqrt(3.0) / 6.0) * std::atan((x - 1.0) / std::sqrt(3.0));
    };

    const auto x_0 = -1.0;
    const auto x_1 = 1.0;
    const auto h_1 = 0.5;
    const auto h_2 = 0.25;

    const double exact_value = exact_integral(x_1) - exact_integral(x_0);

    const std::vector<MethodData> methods = {
        {
            .name = "Left Rectangles",
            .p = 1,
            .func = RectangleLeftMethod
        },
        {
            .name = "Right Rectangles",
            .p = 1,
            .func = RectangleRightMethod
        },
        {
            .name = "Middle Rectangles",
            .p = 2,
            .func = RectangleMiddleMethod
        },
        {
            .name = "Trapezoid",
            .p = 2,
            .func = TrapezoidMethod
        },
        {
            .name = "Simpson",
            .p = 4,
            .func = SimpsonMethod
        }
    };

    std::cout << std::left
              << std::setw(25) << "Method" << std::setw(15) << "I(h_1)"
              << std::setw(15) << "I(h_2)" << std::setw(10) << "Order"
              << std::setw(15) << "Error" << std::setw(15) << "Refined I"
              << std::setw(20) << "Absolute Error" << "\n";
    std::cout << std::string(109, '-') << "\n";

    std::cout << std::left
              << std::setw(25) << "Analytical" << std::setw(15) << exact_value
              << std::setw(15) << exact_value << std::setw(10) << 0
              << std::setw(15) << 0 << std::setw(15) << exact_value
              << std::setw(20) << 0 << "\n";
    for (const auto &method : methods) {
        const auto I_h_1 = method.func(function, x_0, x_1, h_1);
        const auto I_h_2 = method.func(function, x_0, x_1, h_2);
        const auto p = method.p;
        const auto error  = RungeRombergError(I_h_1, I_h_2, method.p);
        const auto I_refined = RungeRombergRefined(I_h_1, I_h_2, method.p);
        const auto abs_error = std::abs(exact_value - I_h_2);

        std::cout << std::left
                  << std::setw(25) << method.name << std::setw(15) << I_h_1
                  << std::setw(15) << I_h_2 << std::setw(10) << p
                  << std::setw(15) << error << std::setw(15) << I_refined
                  << std::setw(20) << abs_error << "\n";
    }

    return 0;
}