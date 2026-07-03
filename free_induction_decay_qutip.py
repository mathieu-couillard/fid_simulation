import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.integrate import cumulative_trapezoid
from qutip import *

# Global font settings (Unchanged)
plt.rcParams.update(
    {
        "font.size": 10,
        "axes.labelsize": 9,
        "figure.labelsize": 9,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 7,
        "legend.title_fontsize": 9,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial"],
        "mathtext.fontset": "custom",
        "mathtext.rm": "Arial",
        "mathtext.it": "Arial:italic",
        "mathtext.bf": "Arial:bold",
    }
)


# --- 2. Simulation Setup ---
curr_to_delta_b = 353 / 28.043

g_col = 2 * np.pi * 4.5
k_ext = 2 * np.pi * 3.2
k_int = k_ext / 100
kappa_tot = k_ext + k_int
t_drive = 0.132
t_max = 0.5
times = np.linspace(0, t_max, 126 * 2)

N_max = 12
a = tensor(destroy(N_max), qeye(N_max))
b = tensor(qeye(N_max), destroy(N_max))
H_drive_op = 1j * np.sqrt(k_ext) * (a.dag() - a)


def drive_shape(t, args):
    t_ramp_up = 0.024
    t_ramp_down = 0.016
    if t < t_ramp_up:
        return t / t_ramp_up
    elif t_ramp_up <= t <= t_drive:
        return 1.0
    elif t_drive < t <= (t_drive + t_ramp_down):
        return 1.0 - (t - t_drive) / t_ramp_down
    else:
        return 0.0


def run_simulation(gamma_val, g_val=g_col):
    H0 = g_val * (a.dag() * b + a * b.dag())
    H = [H0, [H_drive_op, drive_shape]]
    c_ops = [np.sqrt(kappa_tot) * a, np.sqrt(gamma_val) * b]
    res = mesolve(H, tensor(basis(N_max, 0), basis(N_max, 0)), times, c_ops, [a])
    v_in = np.array([drive_shape(t, None) for t in times])
    v_out = v_in - np.sqrt(k_ext) * res.expect[0]
    return np.abs(v_out)


bare_cavity = run_simulation(0, g_val=0)
norm_ref_peak = np.max(bare_cavity)
norm_ref_area = cumulative_trapezoid(bare_cavity**2, times, initial=0)[-1]

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.37, 4))

fig.subplots_adjust(top=0.985, bottom=0.085, left=0.155, right=0.96)

regime_a = 2 * np.pi * np.linspace(4.5, 12, 6)
regime_b = 2 * np.pi * np.linspace(13, 72, 6)
regime_c = 2 * np.pi * np.linspace(4.5, 72, 21)

for ax, regime, label in zip([ax1, ax2], [regime_a, regime_b], ["a)", "b)"]):
    ax.plot(
        times * 1000,
        (bare_cavity / norm_ref_peak),
        "--",
        color="black",
        linewidth=0.8,
        label="Cavity",
    )
    for gamma in regime:
        trace = run_simulation(gamma)
        ax.plot(
            times * 1000,
            (trace / norm_ref_peak),
            linewidth=0.8,
            label=f"{gamma / (2 * np.pi):.1f} MHz",
        )

    ax.set_xlim(0, 500)
    ax.set_ylim(0, 1.0)
    ax.grid(True)
    ax.set_ylabel("Normalized amplitude")
    ax.axvspan(0, t_drive * 1000, color="lightgray", alpha=0.6)
    ax.text(
        0.03,
        0.92,
        label,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="top",
    )
    ax.legend(loc="upper right", frameon=True)

# SPECIFIC EDIT: Remove x-axis numbers and label for panel (a)
ax1.tick_params(labelbottom=False)
ax1.set_xlabel("")

# Keep x-label for panel (b)
ax2.set_xlabel("Time (ns)", labelpad=1)


plt.show()
