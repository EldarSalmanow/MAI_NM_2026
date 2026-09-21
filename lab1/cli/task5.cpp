#include <lab1/matrix.h>

#include <complex>
#include <fstream>
#include <iostream>


auto QRFactorization(
    const Matrix &A,
    const double epsilon
) -> std::tuple<Matrix, Matrix> {
    using size_type = Matrix::size_type;

    const size_type n = A.Rows();

    Matrix Q = Matrix::Identity(n);
    Matrix R = A;

    for (size_type j = 0; j < n - 1; ++j) {
        Vector v(n, 0.0);

        double norm_x = 0.0;
        for (size_type i = j; i < n; ++i) {
            norm_x += R[i][j] * R[i][j];
        }
        norm_x = std::sqrt(norm_x);

        if (norm_x < epsilon) {
            continue;
        }

        for (size_type i = j; i < n; ++i) {
            v[i] = R[i][j];
        }

        const double sign = (v[j] >= 0.0) ? 1.0 : -1.0;
        v[j] += sign * norm_x;

        double v_norm_sq = 0.0;
        for (size_type i = j; i < n; ++i) {
            v_norm_sq += v[i] * v[i];
        }

        if (v_norm_sq < epsilon) {
            continue;
        }

        for (size_type c = j; c < n; ++c) {
            double dot = 0.0;
            for (size_type r = j; r < n; ++r) {
                dot += v[r] * R[r][c];
            }

            const double scale = 2.0 * dot / v_norm_sq;
            for (size_type r = j; r < n; ++r) {
                R[r][c] -= scale * v[r];
            }
        }

        for (size_type r = 0; r < n; ++r) {
            double dot = 0.0;
            for (size_type c = j; c < n; ++c) {
                dot += Q[r][c] * v[c];
            }

            const double scale = 2.0 * dot / v_norm_sq;
            for (size_type c = j; c < n; ++c) {
                Q[r][c] -= scale * v[c];
            }
        }
    }

    return {Q, R};
}

auto ExtractEigenValues(
    const Matrix &A,
    const double epsilon
) -> VectorCF64 {
    using size_type = Matrix::size_type;

    const size_type n = A.Rows();

    std::vector<std::complex<double>> eigen_values;

    for (size_type i = 0; i < n; ++i) {
        if (i < n - 1 && std::abs(A[i + 1][i]) > epsilon) {
            const double a = A[i][i], b = A[i][i+1], c = A[i+1][i], d = A[i+1][i+1];

            const double trace = a + d;
            const double det = a * d - b * c;

            if (const double discriminant = trace * trace - 4.0 * det; discriminant >= 0) {
                eigen_values.emplace_back((trace + std::sqrt(discriminant)) / 2.0, 0.0);
                eigen_values.emplace_back((trace - std::sqrt(discriminant)) / 2.0, 0.0);
            } else {
                eigen_values.emplace_back(trace / 2.0, std::sqrt(-discriminant) / 2.0);
                eigen_values.emplace_back(trace / 2.0, -std::sqrt(-discriminant) / 2.0);
            }

            ++i;
        } else {
            eigen_values.emplace_back(A[i][i], 0.0);
        }
    }

    return VectorCF64(eigen_values);
}

auto QRAlgorithm(
    const Matrix &A,
    const double epsilon,
    const std::uint64_t max_iterations
) -> VectorCF64 {
    using size_type = VectorCF64::size_type;

    Matrix A_k = A;

    auto prev_eigen_values = ExtractEigenValues(A_k, epsilon);

    for (std::uint64_t iteration = 0; iteration < max_iterations; ++iteration) {
        auto [Q, R] = QRFactorization(A_k, epsilon);

        A_k = R * Q;

        auto current_eigen_values = ExtractEigenValues(A_k, epsilon);

        if (current_eigen_values.Size() != prev_eigen_values.Size()) {
            prev_eigen_values = current_eigen_values;

            continue;
        }

        bool converged = true;
        for (size_type i = 0; i < current_eigen_values.Size(); ++i) {
            if (std::abs(current_eigen_values[i] - prev_eigen_values[i]) > epsilon) {
                converged = false;

                break;
            }
        }

        if (converged) {
            return current_eigen_values;
        }

        prev_eigen_values = current_eigen_values;
    }

    return prev_eigen_values;
}


int main(int argc, char **argv) {
    // 1. Исходную матрицу A.
    // 2. Матрицы Q и R.
    // 3. Результат умножения Q · R.
    // 4. Найденные собственные значения.

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

    // 1. print A
    std::cout << "--- Matrix A ---\n" << A << std::endl;

    auto [Q, R] = QRFactorization(A, epsilon);

    // 2. print Q and R
    std::cout << "--- Matrix Q and R ---\n" << Q << "\n\n" << R << std::endl;

    // 3. print Q * R
    std::cout << "--- Matrix Q * R ---\n" << Q * R << std::endl;

    const auto Lambda = QRAlgorithm(A, epsilon, max_iterations);

    // 4. print Lambda
    std::cout << "--- Eigen values ---\n" << Lambda << std::endl;

    return 0;
}
