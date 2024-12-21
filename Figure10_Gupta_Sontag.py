
# Figure 10 in the paper "Cumulative dose responses for adapting biological systems" by Ankit Gupta and Eduardo Sontag


import numpy as np
from scipy.integrate import solve_ivp, cumulative_trapezoid, quad
import matplotlib.pyplot as plt

save_fig = True

# specify parameter x0
x0 = 1.0


def ode_system(t, state, p):
    log_z, y, int_y = state
    z = np.exp(log_z)
    dlog_zdt = - y
    dydt = z - p - y
    dint_ydt = y
    return [dlog_zdt, dydt, dint_ydt]


T_list = [25, 30]
p_list = [1.0, 1.2, 1.5, 2.0, 10.0]
u_r = [550, 800, 1200, 2200, 50000, 200000]
u_l = [400, 550, 600, 1000, 10000, 50000]
fig, axes = plt.subplots(nrows=len(p_list), ncols=len(T_list), figsize=(12, 12))
if axes.ndim == 1:
    axes = axes[:, np.newaxis]

for p in p_list:
    for T in T_list:
        t_span = (0, T)  # time span for the simulation
        rho = (T / max(T_list)) ** 2
        u_range = np.linspace(u_l[p_list.index(p)] * rho, u_r[p_list.index(p)] * rho,
                              min(int(u_r[p_list.index(p)] * rho - u_l[p_list.index(p)] * rho), 10000))
        integral_vals = np.zeros_like(u_range)
        print(f"p = {p}, T = {T} Completed!")
        for u in u_range:
            initial_condition = [u / x0, p, 0.0]
            sol = solve_ivp(ode_system, t_span, [np.log(u / x0), 0, 0], args=(p,), method='RK45', rtol=1e-10,
                            atol=1e-12)
            integral_vals[u_range == u] = sol.y[2][-1]
        axes[p_list.index(p), T_list.index(T)].plot(u_range, integral_vals)
        axes[p_list.index(p), T_list.index(T)].set_title(f"p = {p}, T = {T}")
        axes[p_list.index(p), T_list.index(T)].grid(True)

fig.supxlabel('u', fontsize=20)
fig.supylabel('$\int_0^T \\tilde{y}(t)dt$', fontsize=20)
# Adjust layout for better spacing
plt.tight_layout()
if save_fig:
    plt.savefig('NFB_Table.pdf')
# Show the plots
plt.show()
