#include "casadi_mpc/casadi_mpc_template.hpp"
#include "casadi_mpc/matplotlibcpp.h"
#include <casadi/casadi.hpp>
#include <iostream>
#include <string>

#include <chrono>
#include <ros/ros.h>
class CartpoleProb : public casadi_mpc_template::Problem
{
  public:
    CartpoleProb() : Problem(DynamicsType::ContinuesRK4, 4, 1, 60, 0.001)
    {
        using namespace casadi;
        x_ref = {0, M_PI, 0, 0};
        Q = DM::diag({5.0, 30.0, 0.01, 0.01});
        R = DM::diag({0.1});
        Qf = DM::diag({100.0, 100.0, 0.01, 0.01});

        // ub is problem
        Eigen::VectorXd u_lb = (Eigen::VectorXd(1) << -2.0).finished();
        Eigen::VectorXd u_Ub = (Eigen::VectorXd(1) << 2.0).finished();
        Eigen::VectorXd x_lb = (Eigen::VectorXd(4) << -2.0, -2 * M_PI, -50, -50).finished();
        Eigen::VectorXd x_Ub = (Eigen::VectorXd(4) << 2.0, 2 * M_PI, 50, 50).finished();

        set_input_bound(u_lb, u_Ub);
        set_state_bound(x_lb, x_Ub);
    }

    virtual casadi::MX dynamics(casadi::MX x, casadi::MX u) override
    {
        using namespace casadi;
        auto x_cart = x(0);
        auto th = x(1);
        auto x_cart_dot = x(2);
        auto dth = x(3);

        auto d_x_cart_dot =
            (-m * l * sin(th) * dth * dth + m * g * cos(th) * sin(th) + u) / (mcart + m - m * cos(th) * cos(th));
        auto d_dth = (-m * l * cos(th) * sin(th) * dth * dth + u * cos(th) + (mcart + m) * g * sin(th)) /
                     (l * (mcart + m - m * cos(th) * cos(th)));

        return vertcat(x_cart_dot, dth, d_x_cart_dot, d_dth);
    }

    Eigen::VectorXd discretized_dynamics(double dt, Eigen::VectorXd x, Eigen::VectorXd u)
    {
        auto dynamics = [&](Eigen::VectorXd x, Eigen::VectorXd u) -> Eigen::VectorXd {
            auto x_cart = x(0);
            auto th = x(1);
            auto x_cart_dot = x(2);
            auto dth = x(3);

            auto d_x_cart_dot =
                (-m * l * sin(th) * dth * dth + m * g * cos(th) * sin(th) + u(0)) / (mcart + m - m * cos(th) * cos(th));
            auto d_dth = (-m * l * cos(th) * sin(th) * dth * dth + u(0) * cos(th) + (mcart + m) * g * sin(th)) /
                         (l * (mcart + m - m * cos(th) * cos(th)));

            return (Eigen::VectorXd(4) << x_cart_dot, dth, d_x_cart_dot, d_dth).finished();
        };

        return casadi_mpc_template::integrate_dynamics_rk4<Eigen::VectorXd>(dt, x, u, dynamics);
    }

    virtual casadi::MX stage_cost(casadi::MX x, casadi::MX u) override
    {
        using namespace casadi;
        MX L = 0;
        auto e = x - x_ref;
        L += 0.5 * mtimes(e.T(), mtimes(Q, e));
        L += 0.5 * mtimes(u.T(), mtimes(R, u));

        return dt() * L;
    }

    virtual casadi::MX terminal_cost(casadi::MX x)
    {
        using namespace casadi;
        auto e = x - x_ref;
        return 0.5 * mtimes(e.T(), mtimes(Qf, e));
    }

    const double l = 2.0;     // [m]
    const double m = 1.0;     // [kg]
    const double g = 9.81;    // [m/s^2]
    const double mcart = 0.5; // [kg]

    casadi::DM x_ref;
    casadi::DM Q, R, Qf;
};

void animate(const std::vector<double> &x, const std::vector<double> &angle, const std::vector<double> &u,
             size_t skip = 4)
{
    namespace plt = matplotlibcpp;
    size_t cnt = 1;
    for (size_t i = 0; i < x.size(); i += skip)
    {
        plt::clf();

        double cart_width = 0.5;
        double cart_height = 0.2;
        double pole_length = 0.5;
        double pole_width = 0.05;

        double cart_x = x[i];
        double cart_y = 0.0;
        double pole_x = cart_x + pole_length * sin(angle[i]);
        double pole_y = cart_y - pole_length * cos(angle[i]);

        std::vector<double> cart_x_data = {cart_x - cart_width / 2, cart_x + cart_width / 2, cart_x + cart_width / 2,
                                           cart_x - cart_width / 2, cart_x - cart_width / 2};
        std::vector<double> cart_y_data = {cart_y - cart_height / 2, cart_y - cart_height / 2, cart_y + cart_height / 2,
                                           cart_y + cart_height / 2, cart_y - cart_height / 2};
        std::vector<double> pole_x_data = {cart_x, pole_x};
        std::vector<double> pole_y_data = {cart_y, pole_y};

        plt::set_aspect(1.0);
        plt::plot(cart_x_data, cart_y_data, "b-");
        plt::plot(pole_x_data, pole_y_data, "r-");
        double x_max = *std::max_element(x.begin(), x.end());
        double x_min = *std::min_element(x.begin(), x.end());
        double range_max = std::max(std::abs(x_max), std::abs(x_min)) + cart_width;
        plt::xlim(-range_max, range_max);
        plt::ylim(-range_max, range_max);
        // plt::save("./movie/sample" + std::to_string(cnt++) + ".png");
        plt::pause(0.001);
        // std::cout << i+1 << "/" << x.size() << std::endl;
    }
}
int main(int argc, char **argv)
{
    // Initialize ROS node
    ros::init(argc, argv, "mpc_cartpole_node");
    ros::NodeHandle nh;

    using namespace casadi_mpc_template;
    std::cout << "MPC Cartpole example" << std::endl;
    auto prob = std::make_shared<CartpoleProb>();
    MPC mpc(prob);

    Eigen::VectorXd x = Eigen::VectorXd::Zero(prob->nx());

    x << 0.3, M_PI, 0, 0; // initial state

    const double dt = 0.001;
    const size_t sim_len = 4000;
    std::vector<double> i_log(sim_len), t_log(sim_len), x_log(sim_len), angle_log(sim_len), u_log(sim_len),
        dt_log(sim_len);

    auto t_all_start = std::chrono::system_clock::now();
    for (size_t i = 0; i < sim_len; i++)
    {
        auto t_start = std::chrono::system_clock::now();

        // Solve for optimal input using MPC
        Eigen::VectorXd u = mpc.solve(x);

        auto t_end = std::chrono::system_clock::now();

        x = prob->discretized_dynamics(dt, x, u);

        double solve_time = std::chrono::duration_cast<std::chrono::microseconds>(t_end - t_start).count() * 1e-6;
        std::cout << "Solve time: " << solve_time << std::endl;
        t_log[i] = i * dt;
        i_log[i] = i;
        x_log[i] = x[0];
        angle_log[i] = x[1];
        u_log[i] = u[0];
        dt_log[i] = solve_time * 1e3;
        std::cout << "state: " << std::endl << x.transpose() << std::endl;
        std::cout << "input: " << std::endl << u.transpose() << std::endl;
    }
    auto t_all_end = std::chrono::system_clock::now();
    double all_time = std::chrono::duration_cast<std::chrono::microseconds>(t_all_end - t_all_start).count() * 1e-6;
    std::cout << "All time: " << all_time << std::endl;

    namespace plt = matplotlibcpp;
    plt::figure();
    plt::named_plot("u", t_log, u_log);
    plt::named_plot("x", t_log, x_log);
    plt::named_plot("angle", t_log, angle_log);
    plt::show();

    plt::figure();
    plt::plot(i_log, dt_log);
    plt::xlabel("iteration");
    plt::ylabel("MPC solve time [ms]");
    plt::show();

    animate(x_log, angle_log, u_log);

    ros::shutdown();

    return 0;
}
