# Cavity-Spin Ensemble Time-Domain Simulation

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/)
[![Framework: QuTiP](https://img.shields.io/badge/Framework-QuTiP-orange.svg)](https://qutip.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)

This repository contains a core time-domain simulation script that models the coherent and dissipative dynamics of a spin ensemble coupled to a single-mode resonator. By solving the master equation under a pulsed drive, this script maps out cavity reflection/transmission amplitudes across different spin dissipation regimes to investigate collective impedance matching conditions.

---

## Physical Model & Implementation

The script utilizes **QuTiP (Quantum Toolbox in Python)** to implement a two-mode Jaynes-Cummings/Tavis-Cummings master equation solver:

- **Hamiltonian:** Includes the standard cavity-spin exchange interaction $H_0 = g_{\text{col}}(a^\dagger b + a b^\dagger)$ combined with a time-dependent, trapezoidal pulse shape driving the cavity mode.
- **Dissipation:** Incorporates cavity decay ($\kappa_{\text{tot}} = \kappa_{\text{ext}} + \kappa_{\text{int}}$) and inhomogeneous spin broadening modeled as an effective spin decay rate ($\gamma$).
- **Regimes Explored:** Sweeps $\gamma$ across multiple parameter sets to demonstrate the transition from undercoupled/overcoupled behavior to the critical impedance matching point.

