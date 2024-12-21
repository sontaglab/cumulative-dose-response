# Figure 8 in the paper "Cumulative dose responses for adapting biological systems" by Ankit Gupta and Eduardo Sontag

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt


x0 = 1.0
y_eq = 1.0  # also the steady state for y. It is the same as p.
T_upper_bound = 100.0
u = 300.0
t_span = (0, T_upper_bound)


def compute_Hamiltonians(y):
    log_z, y1, beta, beta_prime, score = y
    z = np.exp(log_z)
    H_beta = beta ** 2 + (beta_prime ** 2) / z
    H = y1 ** 2 + ((z - y_eq - y1) ** 2) / z
    return H, H_beta


def stop_condition(t, state):
    H, H_beta = compute_Hamiltonians(state)
    return max(H - 4, 0) + max(H_beta - 1, 0)


stop_condition.terminal = False


def ode_system(t, state):
    log_z, y, beta, beta_prime, score = state
    z = np.exp(log_z)
    dlog_zdt = - y
    dydt = z - y_eq - y
    dbeta_dt = beta_prime
    dbeta_prime_dt = - beta_prime - z * beta
    dscore_dt = max(beta - 1, 0)
    return [dlog_zdt, dydt, dbeta_dt, dbeta_prime_dt, dscore_dt]



sol = solve_ivp(ode_system, t_span, [np.log(u / x0), 0, 1, 0, 0], events=stop_condition,method='RK45', rtol=1e-10, atol=1e-12)
T_u = T_upper_bound
if sol.t_events[0].size > 0:
    T_u = sol.t_events[0][0]

t_eval = sol.t
beta_values = sol.y[2]
plt.xlabel("t")
plt.ylabel('$u \int_0^t \partial_u y(s)ds$')
plt.plot(t_eval, 1-beta_values)
plt.axhline(y=0, color='g', linestyle='--')
plt.axvline(x=T_u, color='r', linestyle='--',linewidth=1, label='$T_u$')
#plt.grid(True)
plt.title(f'$x_0$ = {x0}, $p$ = {y_eq}, u = {u}')
plt.legend()
plt.savefig(f'NFB_plot_{u}.pdf')
plt.show()
