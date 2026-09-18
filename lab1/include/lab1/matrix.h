#ifndef MAI_NM_2026_MATRIX_H
#define MAI_NM_2026_MATRIX_H

#include <lab1/vector.h>


template<typename T>
class BaseMatrix {
public:
    using value_type = T;

    using container_type = std::vector<BaseVector<value_type>>;

    using iterator = container_type::iterator;

    using const_iterator = container_type::const_iterator;

    using size_type = container_type::size_type;

public:
    BaseMatrix() = default;

    explicit BaseMatrix(container_type data)
            : data_(std::move(data)) {}

    explicit BaseMatrix(size_type dimension, const value_type &default_value = value_type{})
            : BaseMatrix(dimension, dimension, std::move(default_value)) {}

    explicit BaseMatrix(size_type rows, size_type cols, const value_type &default_value = value_type{})
            : data_(rows, BaseVector<value_type>(cols, default_value)) {}

    BaseMatrix(std::initializer_list<std::initializer_list<value_type>> init)
            : data_(init) {}

public:
    static auto New(container_type data) -> BaseMatrix {
        return BaseMatrix(std::move(data));
    }

    static auto Zero(size_type rows, size_type cols) -> BaseMatrix {
        return BaseMatrix(rows, cols, value_type {});
    }

    static auto Zero(size_type n) -> BaseMatrix {
        return BaseMatrix(n, n, value_type {});
    }

    static auto Diagonal(size_type n, const value_type &value) -> BaseMatrix {
        auto matrix = Zero(n);

        for (size_type i = 0; i < n; ++i) {
            matrix[i][i] = value;
        }

        return matrix;
    }

    static auto Identity(size_type n) -> BaseMatrix {
        return Diagonal(n, value_type { 1 });
    }

public:
    auto SwapRows(size_type row_1, size_type row_2) -> void {
        std::swap(data_[row_1], data_[row_2]);
    }

    auto Transpose() const -> BaseMatrix {
        auto result = BaseMatrix(Cols(), Rows());

        for (size_type i = 0; i < Rows(); ++i) {
            for (size_type j = 0; j < Cols(); ++j) {
                result[j][i] = data_[i][j];
            }
        }

        return result;
    }

    auto NormInfinity() const -> value_type {
        value_type max_sum {};

        for (size_type i = 0; i < Rows(); ++i) {
            value_type current_sum {};

            for (size_type j = 0; j < Cols(); ++j) {
                current_sum += std::abs(data_[i][j]);
            }

            if (current_sum > max_sum) {
                max_sum = current_sum;
            }
        }

        return max_sum;
    }

    auto Rows() const -> size_type {
        return data_.size();
    }

    auto Cols() const -> size_type {
        return Empty() ? 0 : data_[0].Size();
    }

    auto Empty() const -> bool {
        return data_.empty();
    }

public:
    auto operator+() const -> BaseMatrix {
        auto matrix = *this;

        for (size_type i = 0; i < Rows(); ++i) {
            matrix[i] = +matrix[i];
        }

        return matrix;
    }

    auto operator-() const -> BaseMatrix {
        auto matrix = *this;

        for (size_type i = 0; i < Rows(); ++i) {
            matrix[i] = -matrix[i];
        }

        return matrix;
    }

    auto operator+(const BaseMatrix &other) const -> BaseMatrix {
        auto matrix = *this;

        matrix += other;

        return matrix;
    }

    auto operator-(const BaseMatrix &other) const -> BaseMatrix {
        auto matrix = *this;

        matrix -= other;

        return matrix;
    }

    auto operator*(const BaseMatrix &other) const -> BaseMatrix {
        if (Cols() != other.Rows()) {
            throw std::runtime_error("[ERROR] Can`t do multiplication with mismatched matrix dimensions!");
        }

        BaseMatrix result(Rows(), other.Cols());

        for (size_type i = 0; i < Rows(); ++i) {
            for (size_type j = 0; j < other.Cols(); ++j) {
                for (size_type k = 0; k < Cols(); ++k) {
                    result[i][j] += data_[i][k] * other[k][j];
                }
            }
        }

        return result;
    }

    auto operator*(const BaseVector<value_type> &vector) const -> BaseVector<value_type> {
        if (Cols() != vector.Size()) {
            throw std::runtime_error("[ERROR] Can`t do multiplication with mismatched matrix and vector dimensions!");
        }

        BaseVector<value_type> result(Rows(), value_type{});
        for (size_type i = 0; i < Rows(); ++i) {
            for (size_type j = 0; j < Cols(); ++j) {
                result[i] += data_[i][j] * vector[j];
            }
        }

        return result;
    }

    auto operator+=(const BaseMatrix &other) -> BaseMatrix & {
        if (Rows() != other.Rows()
         || Cols() != other.Cols()) {
            throw std::runtime_error("[ERROR] Can`t do addition with different dimensions of matrix!");
        }

        for (size_type i = 0; i < Rows(); ++i) {
            data_[i] += other[i];
        }

        return *this;
    }

    auto operator-=(const BaseMatrix &other) -> BaseMatrix & {
        if (Rows() != other.Rows()
         || Cols() != other.Cols()) {
            throw std::runtime_error("[ERROR] Can`t do substraction with different dimensions of matrix!");
        }

        for (size_type i = 0; i < Rows(); ++i) {
            data_[i] -= other[i];
        }

        return *this;
    }

    auto operator*=(const BaseMatrix &other) -> BaseMatrix & {
        *this = *this * other;

        return *this;
    }

    auto operator[](size_type index) -> BaseVector<value_type> & {
        if (index >= Rows()) {
            throw std::runtime_error("[ERROR] Index out of range!");
        }

        return data_[index];
    }

    auto operator[](size_type index) const -> const BaseVector<value_type> & {
        if (index >= Rows()) {
            throw std::runtime_error("[ERROR] Index out of range!");
        }

        return data_[index];
    }

    friend auto operator<<(std::ostream &stream, const BaseMatrix &matrix) -> std::ostream & {
        for (size_type i = 0; i < matrix.Rows(); ++i) {
            for (size_type j = 0; j < matrix.Cols(); ++j) {
                stream
                    << std::setw(10) << std::fixed << std::setprecision(4)
                    << static_cast<double>(matrix[i][j]) << " ";
            }

            stream << std::endl;
        }

        return stream;
    }

    friend auto operator>>(std::istream &stream, BaseMatrix &matrix) -> std::istream & {
        size_type rows, columns;
        stream >> rows >> columns;

        matrix.data_.assign(rows, BaseVector<value_type>(columns, value_type {}));

        for (size_type i = 0; i < rows; ++i) {
            for (size_type j = 0; j < columns; ++j) {
                stream >> matrix[i][j];
            }
        }

        return stream;
    }

private:
    std::vector<BaseVector<value_type>> data_;
};

using MatrixF64 = BaseMatrix<double>;

using Matrix = MatrixF64;

#endif //MAI_NM_2026_MATRIX_H
