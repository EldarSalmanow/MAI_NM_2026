#include <fstream>
#include <iostream>

#include <lab1/vector.h>


auto TridiagonalSolve(const Vector &a, const Vector &b, const Vector &c, const Vector &d) -> std::tuple<Vector, Vector, Vector> {
    using size_type = Vector::size_type;

    const size_type n = b.Size();

    Vector P(n, 0.0);
    Vector Q(n, 0.0);

    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];

    for (size_type i = 1; i < n; ++i) {
        double denominator = b[i] + a[i] * P[i - 1];

        if (i < n - 1) {
            P[i] = -c[i] / denominator;
        } else {
            P[i] = 0.0;
        }

        Q[i] = (d[i] - a[i] * Q[i - 1]) / denominator;
    }

    Vector x(n, 0.0);
    x[n - 1] = Q[n - 1];

    for (size_type i = n - 1; i-- > 0; ) {
        x[i] = P[i] * x[i + 1] + Q[i];
    }

    return {P, Q, x};
}

int main(int argc, char **argv) {
    if (argc < 2) {
        std::cerr << "[ERROR] Please specify the input file as a command line argument!" << std::endl;
        return 1;
    }

    std::ifstream file(argv[1]);

    if (!file.is_open()) {
        std::cerr << "[ERROR] Can't open input file!" << std::endl;

        return 1;
    }

    using size_type = Vector::size_type;

    size_type n;
    file >> n;

    Vector a(n, 0.0);
    Vector b(n, 0.0);
    Vector c(n, 0.0);
    Vector d(n, 0.0);

    for (size_type i = 0; i < n; ++i) {
        file >> a[i];
    }

    for (size_type i = 0; i < n; ++i) {
        file >> b[i];
    }

    for (size_type i = 0; i < n; ++i) {
        file >> c[i];
    }

    for (size_type i = 0; i < n; ++i) {
        file >> d[i];
    }

    auto [P, Q, x] = TridiagonalSolve(a, b, c, d);

    // 1. print P and Q
    std::cout << "--- P ---\n" << P << std::endl;
    std::cout << "--- Q ---\n" << Q << std::endl;

    // 2. print x
    std::cout << "--- x ---\n" << x << std::endl;

    return 0;
}
