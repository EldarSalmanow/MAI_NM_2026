#ifndef MAI_NM_2026_VECTOR_H
#define MAI_NM_2026_VECTOR_H

#include <cmath>
#include <iomanip>
#include <vector>


template<typename T>
class BaseVector {
public:
    using value_type = T;

    using container_type = std::vector<value_type>;

    using iterator = container_type::iterator;

    using const_iterator = container_type::const_iterator;

    using size_type = container_type::size_type;

public:
    BaseVector() = default;

    explicit BaseVector(container_type data)
            : data_(std::move(data)) {}

    explicit BaseVector(size_type size, const value_type &default_value = value_type{})
            : data_(size, default_value) {}

    BaseVector(std::initializer_list<value_type> init)
            : data_(init) {}

public:
    static auto New() -> BaseVector {
        return BaseVector();
    }

    static auto New(container_type data) -> BaseVector {
        return BaseVector(std::move(data));
    }

    static auto New(size_type size, const value_type &default_value = value_type{}) -> BaseVector {
        return BaseVector(size, std::move(default_value));
    }

    static auto New(std::initializer_list<value_type> init) -> BaseVector {
        return BaseVector(init);
    }

public:
    auto Dot(const BaseVector &other) const -> value_type {
        value_type result {};

        for (size_type i = 0; i < Size(); ++i) {
            result += data_[i] * other[i];
        }

        return result;
    }

    auto Norm() const -> value_type {
        return std::sqrt(Dot(*this));
    }

    auto NormInfinity() const -> value_type {
        value_type max_value {};

        for (size_type i = 0; i < Size(); ++i) {
            if (std::abs(data_[i]) > max_value) {
                max_value = std::abs(data_[i]);
            }
        }

        return max_value;
    }

    auto Size() const -> size_type {
        return data_.size();
    }

    auto Empty() const -> bool {
        return data_.empty();
    }

public:
    auto begin() -> iterator {
        return data_.begin();
    }

    auto end() -> const_iterator {
        return data_.end();
    }

    auto begin() const -> const_iterator {
        return data_.begin();
    }

    auto end() const -> const_iterator {
        return data_.end();
    }

public:
    auto operator+() const -> BaseVector {
        auto vector = *this;

        for (size_type i = 0; i < Size(); ++i) {
            vector[i] = +vector[i];
        }

        return vector;
    }

    auto operator-() const -> BaseVector {
        auto vector = *this;

        for (size_type i = 0; i < Size(); ++i) {
            vector[i] = -vector[i];
        }

        return vector;
    }

    auto operator+(const BaseVector &other) const -> BaseVector {
        auto vector = *this;

        vector += other;

        return vector;
    }

    auto operator-(const BaseVector &other) const -> BaseVector {
        auto vector = *this;

        vector -= other;

        return vector;
    }

    auto operator*(const value_type &scalar) const -> BaseVector {
        auto vector = *this;

        vector *= scalar;

        return vector;
    }

    auto operator+=(const BaseVector &other) -> BaseVector & {
        if (Size() != other.Size()) {
            throw std::runtime_error("[ERROR] Can`t do addition with different sizes of vectors!");
        }

        for (size_type i = 0; i < Size(); ++i) {
            data_[i] += other[i];
        }

        return *this;
    }

    auto operator-=(const BaseVector &other) -> BaseVector & {
        if (Size() != other.Size()) {
            throw std::runtime_error("[ERROR] Can`t do substraction with different sizes of vectors!");
        }

        for (size_type i = 0; i < Size(); ++i) {
            data_[i] -= other[i];
        }

        return *this;
    }

    auto operator*=(const value_type &scalar) -> BaseVector & {
        for (size_type i = 0; i < Size(); ++i) {
            data_[i] *= scalar;
        }

        return *this;
    }

    auto operator[](size_type index) -> value_type & {
        if (index >= Size()) {
            throw std::runtime_error("[ERROR] Index out of range!");
        }

        return data_[index];
    }

    auto operator[](size_type index) const -> const value_type & {
        if (index >= Size()) {
            throw std::runtime_error("[ERROR] Index out of range!");
        }

        return data_[index];
    }

    friend auto operator<<(std::ostream &stream, const BaseVector &vector) -> std::ostream & {
        for (size_type i = 0; i < vector.Size(); ++i) {
            stream
                << std::setw(10) << std::fixed << std::setprecision(4)
                << static_cast<double>(vector[i])
                << std::endl;
        }

        return stream;
    }

    friend auto operator>>(std::istream &stream, BaseVector &vector) -> std::istream & {
        size_type size;
        stream >> size;

        vector.data_.assign(size, value_type{});

        for (size_type i = 0; i < size; ++i) {
            stream >> vector[i];
        }

        return stream;
    }

private:
    container_type data_;
};

using VectorF64 = BaseVector<double>;

using Vector = VectorF64;

#endif //MAI_NM_2026_VECTOR_H