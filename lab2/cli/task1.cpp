#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>


template<typename PhiT>
auto SimpleIterationsMethod(
    PhiT &&phi,
    const double initial_x,
    const double epsilon,
    const std::uint64_t max_iterations
) -> double {
    std::cout << "Simple Iterations Method" << std::endl;
    std::cout << std::setw(5) << "Iter" << " | "
              << std::setw(12) << "x_k" << " | "
              << std::setw(15) << "|x_k - x_{k-1}|" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    double x_curr = initial_x;
    double x_prev = 0.0;
    double error = 0.0;
    std::uint64_t iteration = 0;

    do {
        x_prev = x_curr;
        x_curr = phi(x_prev);
        error = std::abs(x_curr - x_prev);
        ++iteration;

        std::cout << std::setw(5) << iteration << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr << " | "
                  << std::setw(15) << std::scientific << std::setprecision(5) << error << std::endl;
    } while (error >= epsilon && iteration < max_iterations);

    std::cout << std::endl;

    return x_curr;
}

template<typename FuncT, typename DerivativeT>
auto NewtonMethod(
    FuncT &&function,
    DerivativeT &&derivative,
    const double initial_x,
    const double epsilon,
    const std::uint64_t max_iterations
) -> double {
    std::cout << "Newton's Method" << std::endl;
    std::cout << std::setw(5) << "Iter" << " | "
              << std::setw(12) << "x_k" << " | "
              << std::setw(15) << "|x_k - x_{k-1}|" << std::endl;
    std::cout << std::string(40, '-') << std::endl;

    double x_curr = initial_x;
    double x_prev = 0.0;
    double error = 0.0;
    std::uint64_t iteration = 0;

    do {
        x_prev = x_curr;
        x_curr = x_prev - function(x_prev) / derivative(x_prev);
        error = std::abs(x_curr - x_prev);
        ++iteration;

        std::cout << std::setw(5) << iteration << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr << " | "
                  << std::setw(15) << std::scientific << std::setprecision(5) << error << std::endl;
    } while (error >= epsilon && iteration < max_iterations);

    std::cout << std::endl;

    return x_curr;
}

int main(int argc, char **argv) {
    // 2.1. Реализовать методы простой итерации и Ньютона решения нелинейных уравнений в виде программ,
    // задавая в качестве входных данных точность вычислений.
    // С использованием разработанного программного обеспечения найти положительный корень нелинейного уравнения
    // (начальное приближение определить графически).
    // Проанализировать зависимость погрешности вычислений от количества итераций

    const auto function = [] (const double x) -> double {
        return std::pow(3, x) - 5 * std::pow(x, 2) + 1;
    };
    const auto phi = [] (const double x) -> double {
        return std::sqrt((std::pow(3, x) + 1.0) / 5.0);
    };
    const auto derivative = [] (const double x) -> double {
        return std::log(3) * std::pow(3, x) - 10 * x;
    };

    if (argc < 2) {
        std::cerr << "[ERROR] Please specify the input file as a command line argument!" << std::endl;

        return 1;
    }

    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "[ERROR] Can`t open input file!" << std::endl;

        return 1;
    }

    double initial_x_1, initial_x_2;
    double epsilon;
    std::uint64_t max_iterations;
    file >> initial_x_1 >> initial_x_2 >> epsilon >> max_iterations;

    // 1. print Simple Iterations Method root
    const double root_simple = SimpleIterationsMethod(phi, initial_x_1, epsilon, max_iterations);
    std::cout << "Simple Iterations Method: x = "
              << std::fixed << std::setprecision(7) << root_simple << std::endl << std::endl;

    // 2. print Newton Method root
    const double root_newton = NewtonMethod(function, derivative, initial_x_2, epsilon, max_iterations);
    std::cout << "Newton Method: x = "
              << std::fixed << std::setprecision(7) << root_newton << std::endl << std::endl;

    return 0;
}
