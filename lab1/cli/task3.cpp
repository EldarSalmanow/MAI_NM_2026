#include <fstream>
#include <iostream>

#include <lab1/matrix.h>


auto TransformToEquivalentForm(const Matrix &A, const Vector &b) -> std::tuple<Matrix, Vector> {
    using size_type = Matrix::size_type;

    const size_type n = A.Rows();

    auto alpha = Matrix::Zero(n);
    Vector beta(n, 0.0);

    for (size_type i = 0; i < n; ++i) {
        beta[i] = b[i] / A[i][i];

        for (size_type j = 0; j < n; ++j) {
            if (i != j) {
                alpha[i][j] = -A[i][j] / A[i][i];
            }
        }
    }

    return {alpha, beta};
}

auto SimpleIteration(const Matrix &alpha, const Vector &beta, const double epsilon, const std::uint64_t max_iterations) -> std::tuple<Vector, std::size_t> {
    const double norm_alpha = alpha.NormInfinity();
    const double eps_coefficient = (norm_alpha < 1.0) ? (1.0 - norm_alpha) / norm_alpha : 1.0;

    Vector x_prev = beta;
    Vector x_curr = beta;
    std::size_t iterations = 0;

    do {
        x_prev = x_curr;

        x_curr = (alpha * x_prev) + beta;

        ++iterations;
    } while ((x_curr - x_prev).NormInfinity() > epsilon * eps_coefficient
          && iterations < max_iterations);

    return {x_curr, iterations};
}

auto SeidelMethod(const Matrix &alpha, const Vector &beta, const double epsilon, const std::uint64_t max_iterations) -> std::tuple<Vector, std::size_t> {
    using size_type = Matrix::size_type;

    const size_type n = alpha.Rows();

    auto C = Matrix::Zero(n);
    for (size_type i = 0; i < n; ++i) {
        for (size_type j = i + 1; j < n; ++j) {
            C[i][j] = alpha[i][j];
        }
    }

    const double norm_alpha = alpha.NormInfinity();
    const double norm_C = C.NormInfinity();
    const double eps_coefficient = (norm_alpha < 1.0) ? (1.0 - norm_alpha) / norm_C : 1.0;

    Vector x_prev = beta;
    Vector x_curr = beta;
    std::size_t iterations = 0;

    do {
        x_prev = x_curr;

        for (size_type i = 0; i < n; ++i) {
            double sum = 0.0;

            for (size_type j = 0; j < i; ++j) {
                sum += alpha[i][j] * x_curr[j];
            }

            for (size_type j = i + 1; j < n; ++j) {
                sum += alpha[i][j] * x_prev[j];
            }

            x_curr[i] = sum + beta[i];
        }

        ++iterations;
    } while ((x_curr - x_prev).NormInfinity() > epsilon * eps_coefficient
          && iterations < max_iterations);

    return {x_curr, iterations};
}

int main(int argc, char **argv) {
    const std::uint64_t max_iterations = 1000000;

    if (argc < 2) {
        std::cerr << "[ERROR] Please specify the input file as a command line argument!" << std::endl;

        return 1;
    }

    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "[ERROR] Can`t open input file!" << std::endl;

        return 1;
    }

    Matrix A;
    Vector b;
    double epsilon;
    file >> A >> b >> epsilon;

    auto [alpha, beta] = TransformToEquivalentForm(A, b);

    // 1. print epsilon
    std::cout << "--- Epsilon ---\n" << epsilon << std::endl;

    double norm_alpha = alpha.NormInfinity();
    std::cout << "Norm ||Alpha|| = " << norm_alpha << std::endl;
    if (norm_alpha >= 1.0) {
        std::cout << "[WARNING] Norm >= 1!" << std::endl << std::endl;
    } else {
        std::cout << "[INFO] Norm < 1." << std::endl << std::endl;
    }

    // 2 & 3 & 4. print iterations, x and speed
    auto [x_simple, iterations_simple] = SimpleIteration(alpha, beta, epsilon, max_iterations);
    std::cout << "=== SIMPLE ITERATION METHOD ===\n";
    std::cout << "Iterations: " << iterations_simple << "\n";
    std::cout << "Solution:\n" << x_simple << "\n";

    auto [x_seidel, iterations_seidel] = SeidelMethod(alpha, beta, epsilon, max_iterations);
    std::cout << "=== SEIDEL METHOD ===\n";
    std::cout << "Iterations: " << iterations_seidel << "\n";
    std::cout << "Solution:\n" << x_seidel << "\n";

    return 0;
}
