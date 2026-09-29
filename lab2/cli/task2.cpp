#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>

#include <lab2/matrix.h>


template<typename PhiT>
auto SimpleIterationsMethod(
    PhiT &&phi,
    const Vector &initial_x,
    const double epsilon,
    const std::uint64_t max_iterations
) -> Vector {
    std::cout << "Simple Iterations Method" << std::endl;
    std::cout << std::setw(5) << "Iter" << " | "
              << std::setw(12) << "x_1" << " | "
              << std::setw(12) << "x_2" << " | "
              << std::setw(15) << "||X_k - X_{k-1}||" << std::endl;
    std::cout << std::string(55, '-') << std::endl;

    Vector x_curr = initial_x;
    Vector x_prev = initial_x;
    std::uint64_t iteration = 0;
    double error = 0.0;

    do {
        x_prev = x_curr;
        x_curr = phi(x_prev);
        error = (x_curr - x_prev).NormInfinity();
        ++iteration;

        std::cout << std::setw(5) << iteration << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr[0] << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr[1] << " | "
                  << std::setw(15) << std::scientific << std::setprecision(5) << error << std::endl;
    } while (error >= epsilon && iteration < max_iterations);

    std::cout << std::endl;

    return x_curr;
}

template<typename FuncT, typename JacobianT>
auto NewtonMethod(
    FuncT &&function,
    JacobianT &&jacobian,
    const Vector &initial_x,
    const double epsilon,
    const std::uint64_t max_iterations
) -> Vector {
    std::cout << "Newton's Method" << std::endl;
    std::cout << std::setw(5) << "Iter" << " | "
              << std::setw(12) << "x_1" << " | "
              << std::setw(12) << "x_2" << " | "
              << std::setw(15) << "||X_k - X_{k-1}||" << std::endl;
    std::cout << std::string(55, '-') << std::endl;

    Vector x_curr = initial_x;
    Vector x_prev = initial_x;
    std::uint64_t iteration = 0;
    double error = 0.0;

    do {
        x_prev = x_curr;

        auto J = jacobian(x_prev);
        auto F = function(x_prev);

        const double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];

        auto J_adj = Matrix::New({
            { J[1][1], -J[0][1] },
            {-J[1][0],  J[0][0] }
        });

        auto dX = (J_adj * (-F)) * (1.0 / det);

        x_curr = x_prev + dX;
        error = dX.NormInfinity();
        ++iteration;

        std::cout << std::setw(5) << iteration << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr[0] << " | "
                  << std::setw(12) << std::fixed << std::setprecision(7) << x_curr[1] << " | "
                  << std::setw(15) << std::scientific << std::setprecision(5) << error << std::endl;
    } while (error >= epsilon && iteration < max_iterations);

    std::cout << std::endl;

    return x_curr;
}

int main(int argc, char **argv) {
    // 2.2. Реализовать методы простой итерации и Ньютона решения систем нелинейных уравнений в виде программного кода,
    // задавая в качестве входных данных точность вычислений.
    // С использованием разработанного программного обеспечения решить систему нелинейных уравнений
    // (при наличии нескольких решений найти то из них, в котором значения неизвестных являются положительными);
    // начальное приближение определить графически.
    // Проанализировать зависимость погрешности вычислений от количества итераций.

    if (argc < 2) {
        std::cerr << "[ERROR] Please specify the input file as a command line argument!" << std::endl;

        return 1;
    }

    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "[ERROR] Can`t open input file!" << std::endl;

        return 1;
    }

    double initial_x_1_1, initial_x_2_1;
    double initial_x_1_2, initial_x_2_2;
    double epsilon;
    std::uint64_t max_iterations;
    file >> initial_x_1_1 >> initial_x_2_1 >> initial_x_1_2 >> initial_x_2_2 >> epsilon >> max_iterations;

    auto initial_x_1 = Vector::New({
        initial_x_1_1,
        initial_x_2_1
    });
    auto initial_x_2 = Vector::New({
        initial_x_1_2,
        initial_x_2_2
    });

    const auto F = [] (const Vector &x) -> Vector {
        return Vector::New({
            x[0] - std::cos(x[1]) - 3.0,
            x[1] - std::sin(x[0]) - 3.0
        });
    };

    const auto J = [] (const Vector &x) -> Matrix {
        return Matrix::New({
            {1.0, std::sin(x[1])},
            {-std::cos(x[0]), 1.0}
        });
    };

    const auto Phi = [] (const Vector &x) -> Vector {
        return Vector::New({
            std::cos(x[1]) + 3.0,
            std::sin(x[0]) + 3.0
        });
    };

    // 1. print Simple Iterations Method root
    const Vector root_simple = SimpleIterationsMethod(Phi, initial_x_1, epsilon, max_iterations);
    std::cout << "Simple Iterations Method Roots" << std::endl
              << std::fixed << std::setprecision(7) << root_simple << std::endl;

    // 2. print Newton Method root
    const Vector root_newton = NewtonMethod(F, J, initial_x_2, epsilon, max_iterations);
    std::cout << "Newton Method Roots" << std::endl
              << std::fixed << std::setprecision(7) << root_newton << std::endl;

    return 0;
}
