# Modelling and Visualisation Projects

Welcome to a collection of simulation projects developed as part of my university course in **Modelling and Visualisation**. Each project focuses on simulating different physical systems or phenomena using computational methods in Python.

This repository showcases my ability to:

- Model real-world systems using numerical techniques
- Implement simulations from first principles
- Visualize and analyze results clearly
- Write clean, documented Python code

---

## 📁 Project Overviews

### `cp1_ising_model` – *Ising Model for Ferromagnetism*

Simulates the 2D Ising model using the Metropolis algorithm to study phase transitions and magnetization behavior.

**Highlights:**
- Monte Carlo Markov Chain (MCMC) simulation
- Energy and magnetization plots vs temperature
- Demonstrates phase transition behavior
- Efficient lattice state updates and animations

---

### `cp2_GOL_and_SIRS` – *Game of Life & SIRS Epidemic Model*

Two simulations in one folder:
- **Game of Life:** A classic cellular automaton to simulate emergent complexity from simple rules.
- **SIRS Model:** A stochastic epidemiological model incorporating Susceptible-Infected-Recovered-Susceptible dynamics.

**Highlights:**
- Demonstrates both deterministic (GOL) and stochastic (SIRS) systems
- Explores oscillatory epidemic behavior
- Interactive grid evolution visualizations

---

### `cp3_CahnHilliard_and_Poisson` – *Phase Separation & Electrostatics*

- **Cahn-Hilliard Equation:** Simulates phase separation in binary mixtures.
- **Poisson Equation:** Solves electrostatics problems using numerical PDE techniques.

**Highlights:**
- Numerical solution of partial differential equations (finite difference scheme)
- Physics-based simulation of domain coarsening
- Efficient solvers and clear visualization of potential fields

---






### 📘 `exam` – *Cahn-Hilliard Equation with Advection* (Exam Project, Produced Simulation from Scratch in 3 Hours)

This project extends the **Cahn-Hilliard phase separation model** by incorporating **advection**, simulating the effects of a velocity field on phase dynamics. Developed as part of a final exam, it demonstrates advanced numerical simulation skills and scientific visualization.

#### Features
- Solves the **Cahn-Hilliard equation** using central finite differences and `convolve2d` from `scipy`
- Implements **advection** with a sinusoidal velocity field
- Produces **animated 2D visualizations** of phase evolution
- Computes **variance at steady state** for different initial conditions
- Outputs results to `.txt` for further analysis

#### Simulation Parts
- **Part A**: Baseline Cahn-Hilliard simulation with animation
- **Part B**: Visual inspection and discussion (see included figures and `README.txt`)
- **Part C**: Plot of steady-state variance vs. initial φ₀
- **Part D**: Advection-enabled simulation to study flow impact






##  Tech & Tools

- Python (NumPy, Matplotlib, SciPy)
- Jupyter Notebooks for interactive development
- Custom plotting and animation tools
- Numerical methods: finite differences, Monte Carlo, lattice updates

---

## 🎯 Skills Demonstrated

- Computational physics & simulation design
- Algorithm development & optimization
- Scientific visualization and data analysis
- Clean code structure and documentation

---

