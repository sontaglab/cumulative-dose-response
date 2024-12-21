# Figure 9 in the paper "Cumulative dose responses for adapting biological systems" by Ankit Gupta and Eduardo Sontag

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

x0 = 1.0
y_eq = 1.0  # It is the steady state for y. Also the same as p.
T_upper_bound = 30.0
u_min = 100.0
u_max = 400.0
number_u_values = 1000

# Define the ODE system
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


stop_condition.terminal = True


def ode_system(t, state):
    log_z, y, beta, beta_prime, score = state
    z = np.exp(log_z)
    dlog_zdt = - y
    dydt = z - y_eq - y
    dbeta_dt = beta_prime
    dbeta_prime_dt = - beta_prime - z * beta
    dscore_dt = max(beta - 1, 0)
    return [dlog_zdt, dydt, dbeta_dt, dbeta_prime_dt, dscore_dt]


u_values = np.linspace(u_min, u_max, number_u_values)
T_u_values = np.zeros_like(u_values)
mono_score = np.zeros_like(u_values)
for u in u_values:
    sol = solve_ivp(ode_system, t_span, [np.log(u / x0), 0, 1, 0, 0], events=stop_condition,method='RK45', rtol=1e-10, atol=1e-12)
    T_u = T_upper_bound
    if sol.t_events[0].size > 0:
        T_u = sol.t_events[0][0]

    t_eval = sol.t
    beta_values = sol.y[2]
    T_u_values[u_values == u] = T_u
    mono_score[u_values == u] = sol.y[4][-1]
    print(f"u = {u}, T = {T_u}", f"non_mono_score = {mono_score[u_values == u]}")

plt.figure(figsize=(12, 12))
plt.subplot(2, 1, 1)
plt.plot(u_values, mono_score)
plt.xlabel("u")
plt.ylabel("Non-monotonicity score $S_u$")
plt.grid(True)

plt.subplot(2, 1, 2)
plt.plot(u_values, T_u_values)
plt.xlabel("u")
plt.ylabel("$T_u$")
plt.grid(True)
#save the figure as a pdf
plt.savefig('NFB_score_tu.pdf')
plt.show()
