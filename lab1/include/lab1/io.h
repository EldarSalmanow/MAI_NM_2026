#ifndef MAI_NM_2026_IO_H
#define MAI_NM_2026_IO_H

#include <istream>

#include <lab1/matrix.h>


auto ReadVector(std::istream &stream) -> Vector;

auto ReadMatrix(std::istream &stream) -> Matrix;

#endif //MAI_NM_2026_IO_H