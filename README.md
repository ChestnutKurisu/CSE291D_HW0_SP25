# HW 0: Double and N-Pendulum Simulation (CSE 291D Physics Simulation – SP25)

This repository contains two physics simulations developed in Python using the **VPython** library: 
1. A **double pendulum** (two-link) simulation.  
2. A more general **N-pendulum** (multi-link) simulation.

Both serve as demonstrations of chaotic motion, visualized in real-time with a custom pastel color palette.

## Table of Contents

- [Overview](#overview)
- [Features and Highlights](#features-and-highlights)
- [Project Structure](#project-structure)
- [Installation and Setup](#installation-and-setup)
- [Execution Instructions](#execution-instructions)
- [Simulation Parameters](#simulation-parameters)
- [References and Resources](#references-and-resources)
- [Contact](#contact)

## Overview

A pendulum with multiple links (masses) is a classic example of a physical system that can exhibit highly complex and chaotic behavior. This project implements:
1. A **double pendulum** with explicit equations of motion.  
2. An **N-pendulum** that constructs the mass matrix and forcing terms from summation-based formulas (following [Yesilyurt 2020][yesilyurt2020]) for an arbitrary number of rods and masses.

**Numerical Integration:** Both scripts use a 4th-order Runge–Kutta (RK4) method to solve the coupled ordinary differential equations (ODEs) of motion.

**Real-Time 3D Rendering:** VPython is used to render the system in an interactive 3D environment.

**Trail Visualization:** Shows the path (trail) of the final mass in each simulation to visually emphasize the chaotic nature of the motion.

**Custom Pastel Palette:** Implements a softer color scheme for rods, masses, and background.

These simulations were developed by **Param Somane** (psomane@ucsd.edu) for CSE 291D at UC San Diego in Spring 2025 (SP25). 

## Features and Highlights

- **Runge-Kutta 4th Order Integrator (RK4):** Ensures accurate solutions of the pendulums’ ODEs.
- **High-Fidelity Rendering:** Smooth 3D animations using VPython, with adjustable frame update rates.
- **Customizable Initial Conditions:** Easily modify angles, lengths, and masses to explore chaotic behavior.
- **Trail Retention:** Displays a long trail of the last mass to highlight patterns (or lack thereof).
- **Pastel Color Palette:** Soothing color scheme for rods, masses, and background.

## Project Structure

```
animation/
│   ├── double_pendulum.mp4
│   └── N_pendulum.mp4
tex/
│   ├── report.tex
│   └── UCSD_CSE_291D_Homework_0.pdf
double_pendulum.py
N_pendulum.py
README.md
```

- **double_pendulum.py**: Main Python script for simulating the classic double pendulum.
- **N_pendulum.py**: Python script generalizing to an N-link pendulum system.
- **animation/**: Contains recorded or exported videos (`double_pendulum.mp4`, `N_pendulum.mp4`).
- **tex/**: Contains the LaTeX source (`report.tex`) and compiled PDF for the homework report.
- **README.md**: Project instructions and documentation.

## Installation and Setup

### Prerequisites

- **Python 3.7+** (tested on Python 3.10)
- **pip** (or any alternative Python package manager)
- A functioning **OpenGL**-capable environment is recommended since VPython uses 3D rendering.

## Execution Instructions

1. **Activate the Virtual Environment** (if using one).

2. **Run the Double Pendulum Simulation**:
   ```bash
   python double_pendulum.py
   ```
   A VPython window will open, displaying a 3D animation of the two-link pendulum.

3. **Run the N-Pendulum Simulation**:
   ```bash
   python N_pendulum.py
   ```
   A VPython window will open, displaying the multi-link pendulum with default `N=8`. You can change `N` (and other parameters) directly in the script.

## Simulation Parameters

**Double Pendulum (double_pendulum.py):**
- **Masses (m1, m2):** Defaults are 2.0 kg and 1.0 kg.
- **Lengths (L1, L2):** Both rods default to 1.5 m.
- **Gravitational Acceleration (g):** 9.81 m/s^2.
- **Initial Angles (theta1, theta2):** Set to 1.2 and 2.1 radians, respectively.
- **Time Step (dt):** 0.0005 s, with `N_SUBSTEPS = 10`.
- **Simulation Duration:** 40.0 seconds.
- **Frame Capture:** Optionally set a capture flag if you wish to save frames.

**N-Pendulum (N_pendulum.py):**
- **Number of bobs (N):** Defaults to 8. 
- **Masses (m_i):** All set to 1.0 kg by default, but can be an array of length `N`.
- **Lengths (L_i):** All rods default to 1.5 m, can be customized similarly.
- **Gravitational Acceleration (g):** 9.81 m/s^2.
- **Initial Angles:** Provided as an array `init_thetas`. By default, a spread of angles (like 0.5*(i+1) for each link).
- **Time Step (dt):** 0.001 s, with `N_SUBSTEPS = 10`.
- **Simulation Duration:** 30.0 seconds.

Feel free to experiment by editing the parameters in each script to see how small changes can lead to wildly different pendulum motion.

## References and Resources

1. **Double Pendulum (Mathematical Background):**
   - T. L. Reber, *Chaotic Dynamics of the Double Pendulum*, American Journal of Physics, 2024.
   - R. Fitzpatrick, *Analysis of a Double Pendulum*, University of Texas, 2018.  

2. **N-Pendulum (General Formulation):**  
   - B. Yesilyurt, *Equations of Motion Formulation of a Pendulum Containing N-point Masses*. [arXiv:1910.12610][yesilyurt2020].

3. **VPython Documentation:**  
   - [https://vpython.org](https://vpython.org)

4. **Runge-Kutta Methods:**  
   - E. Hairer, S. P. Nørsett, G. Wanner, *Solving Ordinary Differential Equations I*, Springer, 2nd ed., 2009.

[yesilyurt2020]: https://arxiv.org/abs/1910.12610

## Contact

For questions, feedback, or collaboration inquiries regarding this project:

**Param Somane**  
Email: [psomane@ucsd.edu](mailto:psomane@ucsd.edu)  
**CSE 291D Physics Simulation – Spring 2025**  
University of California, San Diego
