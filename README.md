# Fast SIM: High-Performance Power Electronics Simulator

[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)]()
[![Python](https://img.shields.io/badge/python-3.10%2B-blue)](https://www.python.org/)
[![C++](https://img.shields.io/badge/C%2B%2B-17-blue)](https://isocpp.org/)
[![License](https://img.shields.io/badge/license-MIT-green)]()

**Fast SIM** is a custom-built circuit simulator designed specifically for the rigorous demands of **Power Electronics**.

Unlike general-purpose simulators (like SPICE or PSIM) that can take minutes to settle into steady-state, **Fast SIM** utilizes a hybrid C++/Python architecture and the **Shooting Method** to solve the boundary value problem directly, finding the steady-state in **milliseconds**.

---

## 🚀 Key Features

* **⚡ Hybrid Architecture:** Core numerical engine written in **C++17** for maximum performance, exposed to **Python** via `pybind11` for flexible scripting and analysis.
* **🧮 Shooting Method:** Implements a Newton-Raphson solver in the time domain to find the steady-state solution ($x(T) - x(0) = 0$) instantly, skipping the thermal/electrical startup transient.
* **🚀 Matrix Pre-Computation:** For switched linear circuits, the solver pre-calculates and inverts system matrices for all switch permutations ($2^N$ states) using **OpenMP** parallelization, reducing the time-step cost to simple matrix-vector multiplications.
* **🔍 Stability Analysis:** Native capability to export the **Jacobian** and Monodromy matrices to analyze converter stability (Eigenvalues) and numerical conditioning.
* **📈 Parametric Sweeps:** Python-driven frequency sweeps and sensitivity analysis.

## 🛠️ Technology Stack

* **Core:** C++17, OpenMP (Multi-threading).
* **Linear Algebra:** [Eigen](https://eigen.tuxfamily.org/) (AVX2/FMA optimizations enabled).
* **Bindings:** Pybind11.
* **Build System:** Setuptools & Cibuildwheel (Cross-platform compilation).

## 📦 Installation

### Prerequisites
* Python 3.10+
* C++ Compiler (MSVC on Windows, GCC/Clang on Linux)
* Git

### Building from Source

Since this project relies on high-performance C++ libraries, you must clone recursively and build locally:

```bash
# 1. Clone the repository with submodules (for Eigen)
git clone https://github.com/PA-PB/fast-sim.git
cd fast-sim
git submodule update --init --recursive

# 2. Build and Install (Compiles C++ backend)
pip install . 
```
## Project Structure
```
fast-spice-llc/
├── src/                # C++ High-Performance Backend
│   ├── Solver.cpp      # Newton-Raphson & Matrix Pre-computation
│   ├── Netlist.cpp     # Circuit Graph & Topology Parsing
│   └── bindings.cpp    # Python/C++ Interface
├── fast_spice/         # Python Package
├── extern/             # Third-party libraries (Eigen)
├──Exemples/
│   └── showcase.py         # Full demo script
└── setup.py            # Compilation & Build Configuration
```
## Example 
For example usage consult ```showcase.py```
