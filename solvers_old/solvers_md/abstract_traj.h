#ifndef abstract_traj_H
#define abstract_traj_H

/**
 * @brief in mathematics, traj is the characheristic line, or ODE
 * d_t X_i = f_i({X_k})
 *
 * here we say f is the ode function;
 *
 * some traj need only f function, some also need f', and even need f''. So this
 * file want to provide a standart interface for this ODE.
 */

class ODEfun {
    double;
}


#endif  // abstract_traj_H