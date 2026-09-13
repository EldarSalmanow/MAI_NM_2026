#ifndef MAI_NM_2026_MATRIX_H
#define MAI_NM_2026_MATRIX_H

#include <lab1/vector.h>


template<typename T>
class BaseMatrix {
public:
    using ValueType = T;

public:
    BaseMatrix() = default;

    explicit BaseMatrix(std::vector<BaseVector<T>> data);

public:
    static auto New(std::vector<BaseVector<T>> data) -> BaseMatrix;

    static auto Diagonal(const T &value, std::uint64_t n) -> BaseMatrix;

    static auto One(std::uint64_t n) -> BaseMatrix;

    static auto Zero(std::uint64_t dimension) -> BaseMatrix;

public:
    auto operator+() -> BaseMatrix &;

    auto operator-() -> BaseMatrix &;

    auto operator+(const BaseMatrix &other) -> BaseMatrix &;

    auto operator-(const BaseMatrix &other) -> BaseMatrix &;

    auto operator*(const BaseMatrix &other) -> BaseMatrix &;

    auto operator[](const std::uint64_t &index) -> BaseVector<T> &;

private:
    std::vector<BaseVector<T>> data_;
};

template<typename T>
auto operator+(const BaseMatrix<T> &a, const BaseMatrix<T> &b) -> BaseMatrix<T>;

void f() {
    BaseMatrix<int> a, b;
    a + b;
}

using MatrixI32 = BaseMatrix<std::int32_t>;

using Matrix = MatrixI32;

#endif //MAI_NM_2026_MATRIX_H
