#include <lab4/differential_equations.h>


int main(int argc, char **argv) {
    // 4.1 Реализовать методы Эйлера, Рунге-Кутты и Адамса 4-го порядка в виде программ, задавая в качестве
    // входных данных шаг сетки h. С использованием разработанного программного обеспечения решить задачу Коши
    // для ОДУ 2-го порядка на указанном отрезке. Оценить погрешность численного решения с использованием
    // метода Рунге – Ромберга и путем сравнения с точным решением.

    const auto function = [] (
        double x,
        const Vector &u
    ) -> Vector {
        return Vector{
            u[1],
            (2.0 * x * u[1] - 2.0 * u[0]) / (x * x + 1.0)
        };
    };

    const auto solution = [] (
        double x
    ) -> double {
        return x - x * x + 1;
    };

    const auto a = 0.0 /* = x_0 */, b = 1.0;
    const auto y_0 = 1.0, dy_0 = 1.0;
    const auto h = 0.1, h2 = 0.05;

    auto [x_grid, y_euler, dy_euler] = EulerMethod(function, a, b, y_0, dy_0, h);
    auto [_1, y_rk4, dy_rk4] = RK4Method(function, a, b, y_0, dy_0, h);
    auto [_2, y_adams, dy_adams] = AdamsMethod(function, a, b, y_0, dy_0, h);

    auto [x_grid_2, y_euler_2, _4] = EulerMethod(function, a, b, y_0, dy_0, h2);
    auto [_5, y_rk4_2, _6] = RK4Method(function, a, b, y_0, dy_0, h2);
    auto [_7, y_adams_2, _8] = AdamsMethod(function, a, b, y_0, dy_0, h2);

    std::vector<double> y_exact;
    for (const auto &x : x_grid) {
        y_exact.push_back(solution(x));
    }

    PrintList("x_grid", x_grid);
    PrintList("y_euler", y_euler);
    PrintList("y_rk4", y_rk4);
    PrintList("y_adams", y_adams);
    PrintList("y_exact", y_exact);

    // PrintList("error_euler_exact_h", RungeRombergError(x_grid, y_euler, solution));
    // PrintList("error_euler_exact_h2", RungeRombergError(x_grid_2, y_euler_2, solution));
    // PrintList("error_rk4_exact_h", RungeRombergError(_1, y_rk4, solution));
    // PrintList("error_rk4_exact_h2", RungeRombergError(_5, y_rk4_2, solution));
    // PrintList("error_adams_exact_h", RungeRombergError(_2, y_adams, solution));
    // PrintList("error_adams_exact_h2", RungeRombergError(_7, y_adams_2, solution));

    std::cout << Vector(RungeRombergError(x_grid, y_euler, solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_grid_2, y_euler_2, solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_grid, y_rk4, solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_grid_2, y_rk4_2, solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_grid, y_adams, solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_grid_2, y_adams_2, solution)).NormInfinity() << std::endl;

    auto error_euler_rr = RungeRombergRefined(y_euler, y_euler_2, 1);
    auto error_rk4_rr = RungeRombergRefined(y_rk4, y_rk4_2, 4);
    auto error_adams_rr = RungeRombergRefined(y_adams, y_adams_2, 4);

    std::cout << Vector(error_euler_rr).NormInfinity() << std::endl;
    std::cout << Vector(error_rk4_rr).NormInfinity() << std::endl;
    std::cout << Vector(error_adams_rr).NormInfinity() << std::endl;

    // for (auto &error : error_euler_rr) {
    //     error = std::abs(error);
    // }
    // for (auto &error : error_rk4_rr) {
    //     error = std::abs(error);
    // }
    // for (auto &error : error_adams_rr) {
    //     error = std::abs(error);
    // }
    //
    // PrintList("error_euler_rr_h2", error_euler_rr);
    // PrintList("error_rk4_rr_h2", error_rk4_rr);
    // PrintList("error_adams_rr_h2", error_adams_rr);

    return 0;
}
