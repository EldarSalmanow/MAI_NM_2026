#include <fstream>
#include <iostream>

#include <lab1/matrix.h>


// PLU factorization
auto PLUFactorization(const Matrix &matrix) -> std::tuple<Matrix, Matrix, Matrix, Matrix::size_type> {
    using size_type = Matrix::size_type;

    const size_type n = matrix.Rows();

    Matrix U = matrix;
    Matrix L = Matrix::Identity(n);
    Matrix P = Matrix::Identity(n);
    size_type swaps = 0;

    for (size_type i = 0; i < n; ++i) {
        size_type pivot_row = i;

        double max_value = std::abs(U[i][i]);
        for (size_type k = i + 1; k < n; ++k) {
            if (std::abs(U[k][i]) > max_value) {
                max_value = std::abs(U[k][i]);

                pivot_row = k;
            }
        }

        if (pivot_row != i) {
            U.SwapRows(i, pivot_row);
            P.SwapRows(i, pivot_row);

            for (size_type j = 0; j < i; ++j) {
                std::swap(L[i][j], L[pivot_row][j]);
            }

            ++swaps;
        }

        for (size_type k = i + 1; k < n; ++k) {
            const double factor = U[k][i] / U[i][i];

            L[k][i] = factor;

            for (size_type j = i; j < n; ++j) {
                U[k][j] -= factor * U[i][j];
            }
        }
    }

    return {P, L, U, swaps};
}

// solving SLE by PLU factorization
auto SolveSLE(const Matrix &P, const Matrix &L, const Matrix &U, const Vector &b) -> Vector {
    using size_type  = Matrix::size_type;

    const size_type n = L.Rows();

    auto Pb = P * b;

    Vector y(n, 0.0);
    for (size_type i = 0; i < n; ++i) {
        y[i] = Pb[i];

        for (size_type j = 0; j < i; ++j) {
            y[i] -= L[i][j] * y[j];
        }
    }

    Vector x(n, 0.0);
    for (size_type i = n; i-- > 0;) {
        x[i] = y[i];

        for (size_type j = i + 1; j < n; ++j) {
            x[i] -= U[i][j] * x[j];
        }

        x[i] /= U[i][i];
    }

    return x;
}

// invert matrix by LU factorization
auto InvertMatrix(const Matrix &P, const Matrix &L, const Matrix &U) -> Matrix {
    using size_type = Matrix::size_type;

    const size_type n = L.Rows();

    Matrix A_invert(n, n, 0.0);

    for (size_type i = 0; i < n; ++i) {
        Vector e(n, 0.0);

        e[i] = 1.0;

        Vector column = SolveSLE(P, L, U, e);

        for (size_type j = 0; j < n; ++j) {
            A_invert[j][i] = column[j];
        }
    }

    return A_invert;
}

// determinant of matrix by PLU factorization
auto DeterminantMatrix(const Matrix &U, Matrix::size_type swaps) -> Matrix::value_type {
    using size_type = Matrix::size_type;

    const size_type n = U.Rows();

    double det = 1.0;

    for (size_type i = 0; i < n; ++i) {
        det *= U[i][i];
    }

    return (swaps % 2 != 0) ? -det : det;
}

int main(int argc, char **argv) {
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
    file >> A >> b;

    auto [P, L, U, swaps] = PLUFactorization(A);

    // 1. print L and U
    std::cout << "--- L ---\n" << L << std::endl;
    std::cout << "--- U ---\n" << U << std::endl;

    // 2. print LU
    std::cout << "--- L * U ---\n" << (L * U) << std::endl;

    // 3. print x
    auto x = SolveSLE(P, L, U, b);
    std::cout << "--- x ---\n" << x << std::endl;

    // 4. print A^{-1}
    auto A_invert = InvertMatrix(P, L, U);
    std::cout << "--- A^{-1} ---\n" << A_invert << std::endl;

    // 5. print det(A)
    auto A_det = DeterminantMatrix(U, swaps);
    std::cout << "--- det(A) ---\n" << A_det << std::endl << std::endl;

    // 6. check A * A^{-1} = E
    std::cout << "--- A * A^{-1} ---\n" << (A * A_invert) << std::endl;

    // 7. check PA = LU
    std::cout << "--- PA ---\n" << (P * A) << std::endl;

    return 0;
}
