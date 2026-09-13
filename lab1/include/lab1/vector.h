#ifndef MAI_NM_2026_VECTOR_H
#define MAI_NM_2026_VECTOR_H

#include <cstdint>
#include <vector>


template<typename T>
class BaseVector {
public:
    BaseVector() = default;

private:
    std::vector<T> data_;
};

using VectorI32 = BaseVector<std::int32_t>;

using Vector = VectorI32;

#endif //MAI_NM_2026_VECTOR_H