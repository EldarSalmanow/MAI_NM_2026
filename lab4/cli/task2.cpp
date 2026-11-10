#include <lab4/differential_equations.h>


int main(int argc, char **argv) {
    // 4.2 Реализовать метод стрельбы и конечно-разностный метод решения краевой задачи для ОДУ в виде программ.
    // С использованием разработанного программного обеспечения решить краевую задачу для обыкновенного
    // дифференциального уравнения 2-го порядка на указанном отрезке. Оценить погрешность численного решения
    // с использованием метода Рунге – Ромберга и путем сравнения с точным решением.

    using MathFunc = double(*)(double);

    const MathFunc p = [] (double x) -> double {
        return -1.0 / (x - 1.0);
    };
    const MathFunc q = [] (double x) -> double {
        return 1.0 / (x * (x - 1.0));
    };
    const MathFunc f_x = [] (double x) -> double {
        return 0.0;
    };
    const MathFunc exact_solution = [] (double x) -> double {
        return 2.0 + x + 2.0 * x * std::log(std::abs(x));
    };

    double a = 1.00001, b = 3.0;
    double alpha_0 = 0.0, alpha_1 = 1.0, A = 3.0;
    double beta_0 = 1.0, beta_1 = -3.0, B = -4.0;

    double h = 0.1;
    double h2 = 0.05;

    auto [x_h, y_shooting_h, dy1] = ShootingMethod(p, q, f_x, a, b, alpha_0, alpha_1, A, beta_0, beta_1, B, h);
    auto [x_h2, y_shooting_h2, dy2] = ShootingMethod(p, q, f_x, a, b, alpha_0, alpha_1, A, beta_0, beta_1, B, h2);

    auto [_1, y_fdm_h, dy3] = FiniteDifferenceMethod(p, q, f_x, a, b, alpha_0, alpha_1, A, beta_0, beta_1, B, h);
    auto [_2, y_fdm_h2, dy4] = FiniteDifferenceMethod(p, q, f_x, a, b, alpha_0, alpha_1, A, beta_0, beta_1, B, h2);

    std::vector<double> x_exact, y_exact;

    const auto h_exact = 0.01;
    const auto n_exact = static_cast<std::uint64_t>(std::round((b - a) / h_exact));

    for (std::uint64_t i = 0; i <= n_exact; ++i) {
        const auto x = a + h_exact * static_cast<double>(i);

        x_exact.push_back(x);
        y_exact.push_back(exact_solution(x));
    }

    PrintList("x_h", x_h);
    PrintList("y_sh_h", y_shooting_h);
    PrintList("y_fdm_h", y_fdm_h);

    PrintList("x_h2", x_h2);
    PrintList("y_sh_h2", y_shooting_h2);
    PrintList("y_fdm_h2", y_fdm_h2);

    PrintList("x_exact", x_exact);
    PrintList("y_exact", y_exact);

    // PrintList("error_sh_exact_h", RungeRombergError(x_h, y_shooting_h, exact_solution));
    // PrintList("error_sh_exact_h2", RungeRombergError(x_h2, y_shooting_h2, exact_solution));
    // PrintList("error_fdm_exact_h", RungeRombergError(_1, y_fdm_h, exact_solution));
    // PrintList("error_fdm_exact_h2", RungeRombergError(_2, y_fdm_h2, exact_solution));

    std::cout << Vector(RungeRombergError(x_h, y_shooting_h, exact_solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_h2, y_shooting_h2, exact_solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_h, y_fdm_h, exact_solution)).NormInfinity() << std::endl;
    std::cout << Vector(RungeRombergError(x_h2, y_fdm_h2, exact_solution)).NormInfinity() << std::endl;

    auto error_sh_rr = RungeRombergRefined(y_shooting_h, y_shooting_h2, 4);
    auto error_fdm_rr = RungeRombergRefined(y_fdm_h, y_fdm_h2, 1);

    std::cout << Vector(error_sh_rr).NormInfinity() << std::endl;
    std::cout << Vector(error_fdm_rr).NormInfinity() << std::endl;

    // for (auto &error : error_sh_rr) {
    //     error = std::abs(error);
    // }
    // for (auto &error : error_fdm_rr) {
    //     error = std::abs(error);
    // }
    //
    // PrintList("error_sh_rr_h2", error_sh_rr);
    // PrintList("error_fdm_rr_h2", error_fdm_rr);

    return 0;
}
