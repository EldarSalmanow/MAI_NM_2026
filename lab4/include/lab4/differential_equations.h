#ifndef MAI_NM_2026_DIFFERENTIAL_EQUATIONS_H
#define MAI_NM_2026_DIFFERENTIAL_EQUATIONS_H

#include <lab4/vector.h>

#include <functional>
#include <iomanip>
#include <iostream>


template<typename FunctionT>
auto EulerMethod(
    FunctionT &&function,
    double a /* = x_0 */, double b,
    double y_0, double dy_0,
    double h
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    std::vector<double> x_grid, y_value, dy_value;

    Vector u = {y_0, dy_0};

    for (std::uint64_t i = 0; i <= n; ++i) {
        const auto x = a + h * static_cast<double>(i);

        x_grid.push_back(x);
        y_value.push_back(u[0]);
        dy_value.push_back(u[1]);

        u += h * function(x, u);
    }

    return {x_grid, y_value, dy_value};
}

template<typename FunctionT>
auto RK4Method(
    FunctionT &&function,
    double a /* = x_0 */, double b,
    double y_0, double dy_0,
    double h
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    std::vector<double> x_grid, y_value, dy_value;

    Vector u = {y_0, dy_0};

    for (std::uint64_t i = 0; i <= n; ++i) {
        const auto x = a + h * static_cast<double>(i);

        x_grid.push_back(x);
        y_value.push_back(u[0]);
        dy_value.push_back(u[1]);

        auto k1 = function(x, u);
        auto k2 = function(x + h / 2.0, u + (h / 2.0) * k1);
        auto k3 = function(x + h / 2.0, u + (h / 2.0) * k2);
        auto k4 = function(x + h, u + h * k3);

        u += (h / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
    }

    return {x_grid, y_value, dy_value};
}

template<typename FunctionT>
auto AdamsMethod(
    FunctionT &&function,
    double a /* = x_0 */, double b,
    double y_0, double dy_0,
    double h
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    const auto rk4_result = RK4Method(std::move(function), a, a + 3 * h, y_0, dy_0, h);

    std::vector<double> x_grid = std::get<0>(rk4_result);
    std::vector<double> y_value = std::get<1>(rk4_result);
    std::vector<double> dy_value = std::get<2>(rk4_result);

    std::vector<Vector> function_values;
    for (size_t i = 0; i < 4; ++i) {
        function_values.push_back(function(x_grid[i], {y_value[i], dy_value[i]}));
    }

    Vector u = {y_value.back(), dy_value.back()};

    for (std::uint64_t i = 4; i <= n; ++i) {
        const auto x = a + h * static_cast<double>(i);

        const auto u_next = u + (h / 24.0) * (55.0 * function_values[3] - 59.0 * function_values[2] + 37.0 * function_values[1] - 9.0 * function_values[0]);

        x_grid.push_back(x);
        y_value.push_back(u_next[0]);
        dy_value.push_back(u_next[1]);

        function_values.erase(function_values.begin());
        function_values.push_back(function(x, u_next));

        u = u_next;
    }

    return {x_grid, y_value, dy_value};
}

template<typename FunctionT>
auto ShootingMethod(
    FunctionT &&p, FunctionT &&q, FunctionT &&function,
    double a, double b,
    double alpha_0, double alpha_1, double A,
    double beta_0, double beta_1, double B,
    double h
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    auto system_ODE = [&] (double x, const Vector &u) -> Vector {
        return {u[1], -p(x) * u[1] - q(x) * u[0] + function(x)};
    };

    auto shoot = [&] (double eta) {
        const double y_a  = (alpha_1 != 0.0) ? eta : (A / alpha_0);
        const double dy_a = (alpha_1 != 0.0) ? (A - alpha_0 * eta) / alpha_1 : eta;

        auto [x_grid, y_vector, dy_vector] = RK4Method(system_ODE, a, b, y_a, dy_a, h);

        const double y_b = y_vector.back();
        const double dy_b = dy_vector.back();
        const double phi = beta_0 * y_b + beta_1 * dy_b - B;

        return std::make_tuple(phi, x_grid, y_vector, dy_vector);
    };

    const auto epsilon = 1e-6;

    double eta_0 = 0.0;
    double eta_1 = 1.0;

    auto [phi_0, x_grid_0, y_0, dy_0] = shoot(eta_0);
    auto [phi_1, x_grid_1, y_1, dy_1] = shoot(eta_1);

    std::vector<double> result_x = std::move(x_grid_1);
    std::vector<double> result_y = std::move(y_1);
    std::vector<double> result_dy = std::move(dy_1);

    const auto max_iterations = 100;
    for (int iteration = 0; iteration < max_iterations; ++iteration) {
        if (std::abs(phi_1) <= epsilon) {
            break;
        }

        double eta_2 = eta_1 - ((eta_1 - eta_0) / (phi_1 - phi_0)) * phi_1;

        auto [phi_2, x_grid_2, y_2, dy_2] = shoot(eta_2);

        eta_0 = eta_1;
        phi_0 = phi_1;

        eta_1 = eta_2;
        phi_1 = phi_2;

        result_x = std::move(x_grid_2);
        result_y = std::move(y_2);
        result_dy = std::move(dy_2);
    }

    return {result_x, result_y, result_dy};
}

template<typename FunctionT>
auto FiniteDifferenceMethod(
    FunctionT &&p, FunctionT &&q, FunctionT &&function,
    double a, double b,
    double alpha_0, double alpha_1, double A,
    double beta_0, double beta_1, double B,
    double h
) -> std::tuple<std::vector<double>, std::vector<double>, std::vector<double>> {
    const auto n = static_cast<std::uint64_t>(std::round((b - a) / h));

    std::vector<double> x_grid(n + 1), y_result(n + 1), dy_result(n + 1, 0.0);

    for(int i = 0; i <= n; ++i) {
        x_grid[i] = a + i * h;
    }

    std::vector<double> lower(n + 1), main_d(n + 1), upper(n + 1), rhs(n + 1);

    main_d[0] = alpha_0 - alpha_1 / h;
    upper[0] = alpha_1 / h;
    rhs[0] = A;

    for (int i = 1; i < n; ++i) {
        double x_i = x_grid[i];

        lower[i] = 1.0 - (h / 2.0) * p(x_i);
        main_d[i] = -2.0 + h * h * q(x_i);
        upper[i] = 1.0 + (h / 2.0) * p(x_i);
        rhs[i] = h * h * function(x_i);
    }

    lower[n] = -beta_1 / h;
    main_d[n] = beta_0 + beta_1 / h;
    rhs[n] = B;

    std::vector<double> c(n + 1), d(n + 1);
    c[0] = upper[0] / main_d[0];
    d[0] = rhs[0] / main_d[0];

    for (int i = 1; i <= n; ++i) {
        double denominator = main_d[i] - lower[i] * c[i - 1];

        c[i] = upper[i] / denominator;
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / denominator;
    }

    y_result[n] = d[n];
    for (int i = static_cast<std::int32_t>(n) - 1; i >= 0; --i) {
        y_result[i] = d[i] - c[i] * y_result[i + 1];
    }

    return {x_grid, y_result, dy_result};
}

template<typename FunctionT>
auto RungeRombergError(
    const std::vector<double>& x_grid,
    const std::vector<double>& y_num,
    FunctionT &&function
) -> std::vector<double> {
    std::vector<double> abs_error(x_grid.size(), 0.0);

    for (std::uint64_t i = 0; i < x_grid.size(); ++i) {
        abs_error[i] = std::abs(function(x_grid[i]) - y_num[i]);
    }

    return abs_error;
}

auto RungeRombergRefined(
    const std::vector<double> &y_h,
    const std::vector<double> &y_h2,
    std::uint64_t p
) -> std::vector<double> {
    std::vector<double> error(y_h.size(), 0.0);

    const auto denominator = std::pow(2.0, p) - 1;

    for (std::uint64_t i = 0; i < y_h.size(); ++i) {
        error[i] = (y_h[i] - y_h2[2 * i]) / denominator;
    }

    return error;
}

auto PrintList(
    std::string_view name,
    const std::vector<double> &vector
) -> void {
    std::cout << name << " = [\n";

    for (size_t i = 0; i < vector.size(); ++i) {
        std::cout << std::setprecision(15) << std::fixed << "\t" << vector[i];

        if (i != vector.size() - 1) {
            std::cout << ", ";
        }

        std::cout << "\n";
    }

    std::cout << "]\n";
}

#endif //MAI_NM_2026_DIFFERENTIAL_EQUATIONS_H
