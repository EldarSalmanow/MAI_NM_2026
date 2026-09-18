#include <lab1/matrix.h>

#include <fstream>
#include <iostream>
#include <numbers>


// недиагональная норма Фробениуса
auto OffDiagonalNorm(const Matrix &matrix) -> Matrix::value_type {
    using value_type = Matrix::value_type;
    using size_type = Matrix::size_type;

    value_type result = 0;
    for (size_type i = 1; i < matrix.Rows(); ++i) {
        for (size_type j = 0; j < i; ++j) {
            result += matrix[i][j] * matrix[i][j];
        }
    }

    result = std::sqrt(result);

    return result;
}

auto JacobiAlgorithm(
    Matrix A,
    const double epsilon,
    const std::uint64_t max_iterations
) -> std::tuple<Vector, Matrix> {
    using value_type = Matrix::value_type;
    using size_type = Matrix::size_type;

    auto V = Matrix::Identity(A.Rows());

    std::uint64_t iteration = 1;
    while (OffDiagonalNorm(A) > epsilon && iteration < max_iterations) {
        value_type maximum = 0;
        size_type max_row = 0, max_column = 0;

        for (std::uint64_t i = 1; i < A.Rows(); ++i) {
            for (std::uint64_t j = 0; j < i; ++j) {
                if (const auto element = A[i][j]; std::abs(element) > std::abs(maximum)) {
                    maximum = element;

                    max_row = i;
                    max_column = j;
                }
            }
        }

        auto U = Matrix::Identity(A.Rows());

        double angle = 0.0;

        if (A[max_row][max_row] != A[max_column][max_column]) {
            angle = 1.0 / 2.0 * std::atan(2.0 * maximum / (A[max_row][max_row] - A[max_column][max_column]));
        } else {
            angle = std::numbers::pi / 4.0;
        }

        U[max_row][max_row] = std::cos(angle);
        U[max_row][max_column] = -std::sin(angle);
        U[max_column][max_row] = std::sin(angle);
        U[max_column][max_column] = std::cos(angle);

        A = U.Transpose() * A * U;
        V *= U;

        ++iteration;
    }

    auto Lambda = Vector::New(A.Rows());

    for (size_type i = 0; i < A.Rows(); ++i) {
        Lambda[i] = A[i][i];
    }

    return {Lambda, V};
}


int main(int argc, char **argv) {
    // 1. Заданную точность вычислений ε.
    // 2. Найденные собственные значения.
    // 3. Найденные собственные векторы и матрицу собственных векторов.
    // 4. Проверку A · V = V · Λ, где V — матрица собственных векторов, а Λ — диагональная матрица собственных значений.

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
    double epsilon;
    file >> A >> epsilon;

    // 1. print epsilon
    std::cout << "--- Epsilon ---\n" << epsilon << std::endl;

    auto [Lambda, V] = JacobiAlgorithm(A, epsilon, max_iterations);

    // 2. print eigen values
    std::cout << "--- Eigen values ---\n" << Lambda << std::endl;

    // 3. print eigen vectors
    std::cout << "--- Eigen vectors ---\n" << V << std::endl;

    auto Lambda_matrix = Matrix::Zero(Lambda.Size());
    for (std::uint64_t i = 0; i < Lambda.Size(); ++i) {
        Lambda_matrix[i][i] = Lambda[i];
    }

    // 4. check A * V = V * Lambda
    std::cout << "--- Check A * V = V * Lambda ---\n"
    << "\nA * V\n" << A * V << std::endl
    << "\nV * Lambda\n" << V * Lambda_matrix << std::endl;

    return 0;
}
