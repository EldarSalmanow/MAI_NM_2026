#include <fstream>
#include <iostream>

#include <lab1/io.h>


// PLU factorization
auto PLUFactorization(const Matrix &matrix) -> std::tuple<Matrix, Matrix, Matrix> {
    int n = matrix.size();
    Matrix U = matrix;
    Matrix L = matrix;
    Matrix P = matrix;

    // Инициализация L и P как единичных матриц
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            L[i][j] = (i == j) ? 1.0 : 0.0;
            P[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }

    for (int i = 0; i < n; ++i) {
        // Частичный выбор опорного элемента
        int pivot_row = i;
        double max_val = std::abs(U[i][i]);
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(U[k][i]) > max_val) {
                max_val = std::abs(U[k][i]);
                pivot_row = k;
            }
        }

        // Перестановка строк
        if (pivot_row != i) {
            for(int j = 0; j < n; ++j) std::swap(U[i][j], U[pivot_row][j]);
            for(int j = 0; j < n; ++j) std::swap(P[i][j], P[pivot_row][j]);
            for(int j = 0; j < i; ++j) std::swap(L[i][j], L[pivot_row][j]); // в L меняем только вычисленную часть
        }

        // Исключение по Гауссу
        for (int k = i + 1; k < n; ++k) {
            double factor = U[k][i] / U[i][i];
            L[k][i] = factor;
            for (int j = i; j < n; ++j) {
                U[k][j] -= factor * U[i][j];
            }
        }
    }

    return {P, L, U};
}

// solving SLE by PLU factorization
auto SolveSLE(const Matrix &P, const Matrix &L, const Matrix &U, const Vector &b) -> Vector {

}

// invert matrix by LU factorization
auto InvertMatrix(const Matrix &L, const Matrix &U) -> Matrix;

// determinant of matrix by PLU factorization
auto DeterminantMatrix(const Matrix &P, const Matrix &L, const Matrix &U) -> Matrix::ValueType;

int main(int argc, char **argv) {
    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "[ERROR] Can`t open input file!" << std::endl;

        return 0;
    }

    auto A = ReadMatrix(file);
    auto b = ReadVector(file);

    auto [P, L, U] = PLUFactorization(A);

    // 1. print L and U
    // 2. print LU

    auto x = SolveSLE(P, A, L, U);

    // 3. print x

    auto A_invert = InvertMatrix(L, U);

    // 4. print A^{-1}

    auto A_det = DeterminantMatrix(P, L, U);

    // 5. print det(A)
    // 6. check A * A^{-1} = E
    // 7. check PA = LU

    return 0;
}
