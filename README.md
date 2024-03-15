## CasADi + ROS Model Predictive Control examples
This repository includes examples of implementing MPC using CasADi and ROS.

The Casadi template was based on the repository below.
[simple_casadi_mpc](https://github.com/Kotakku/simple_casadi_mpc/tree/17e26e8cc22eb6400d8714e78163e0fceb23fe9d)

## Dependencies
* Eigen3
* [CasADi](https://github.com/casadi/casadi)
* ROS Noetic
* matplotlib

## Examples
## Inverted Pendulum
```
roscore

rosrun casadi_mpc inverted_pendulum_mpc
```
## Cartpole
```
roscore

rosrun casadi_mpc cartpole_nmpc
```
