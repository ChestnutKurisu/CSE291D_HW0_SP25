# HW 0: Double Pendulum Simulation (CSE 291D Physics Simulation – SP25)

This repository contains a double pendulum physics simulation developed in Python using the **VPython** library. The simulation serves as a demonstration of chaotic motion in a two-link pendulum system, visualized in real-time with a custom pastel color palette.

## Table of Contents

- [Overview](#overview)
- [Features and Highlights](#features-and-highlights)
- [Project Structure](#project-structure)
- [Installation and Setup](#installation-and-setup)
- [Execution Instructions](#execution-instructions)
- [Simulation Parameters](#simulation-parameters)
- [References and Resources](#references-and-resources)
- [Contact](#contact)
- [License](#license)

## Overview

A double pendulum is a classic example of a simple physical system that exhibits highly complex and chaotic behavior. This project implements a **double pendulum simulation** in Python:

1. **Numerical Integration:** Uses a 4th-order Runge–Kutta (RK4) method to solve the coupled ordinary differential equations (ODEs) of motion.
2. **Real-Time 3D Rendering:** VPython is used to render the two-link pendulum in an interactive 3D environment.
3. **Trail Visualization:** Shows the path (trail) of the second mass to visually emphasize the chaotic nature of the motion.
4. **Custom Pastel Palette:** Implements a softer color scheme to create a more visually appealing environment.

This simulation was developed by **Param Somane** (psomane@ucsd.edu) for CSE 291D at UC San Diego in Spring 2025 (SP25). 

## Features and Highlights

- **Runge-Kutta 4th Order Integrator (RK4):** Ensures accurate solutions of the pendulum’s ODEs.
- **High-Fidelity Rendering:** Smooth 3D animations using VPython, with adjustable frame update rates.
- **Customizable Initial Conditions:** Easily modify angles, lengths, and masses to explore chaotic behavior.
- **Trail Retention:** Displays a long trail of the second mass to highlight patterns (or lack thereof).
- **Pastel Color Palette:** Soothing color scheme for rods, masses, and background.

## Project Structure

```
animation/
│   └── double_pendulum.mp4
tex/
│   ├── report.tex
│   └── UCSD_CSE_291D_Homework_0.pdf
double_pendulum.py
README.md
```

- **double_pendulum.py**: Main Python script with the simulation code.
- **animation/**: Contains recorded or exported videos (`double_pendulum.mp4`).
- **tex/**: Contains the LaTeX source and compiled PDF for the homework report.
- **requirements.txt**: Python dependencies (if you are including this file).
- **README.md**: Project instructions and documentation.

## Installation and Setup

### Prerequisites

- **Python 3.7+** (tested on Python 3.10)
- **pip** (or any alternative Python package manager)
- A functioning **OpenGL**-capable environment is recommended since VPython uses 3D rendering.

### Step-by-Step Setup

1. **Clone the Repository** (or download as ZIP):
   ```bash
   git clone https://github.com/ChestnutKurisu/double-pendulum.git
   cd double-pendulum
   ```

2. **Create a Virtual Environment** (recommended):

   Using `venv`:
   ```bash
   python -m venv .venv
   source .venv/bin/activate   # (Linux/macOS)
   .venv\Scripts\activate      # (Windows)
   ```

3. **Install Dependencies**:
   ```bash
   pip install -r requirements.txt
   ```
   *Note: The key libraries are `vpython` and `numpy`. Make sure these install successfully.*

## Execution Instructions

1. **Activate the Virtual Environment** (if using one).
2. **Run the Simulation**:
   ```bash
   python double_pendulum.py
   ```
3. A VPython window will open, displaying the 3D double pendulum. You can rotate/pan/zoom in the window to view the pendulum from different angles.

## Simulation Parameters

- **Masses (m1, m2):** Defaults are 2.0 kg and 1.0 kg.
- **Lengths (L1, L2):** Both rods default to 1.5 m.
- **Gravitational Acceleration (g):** 9.81 m/s^2 (Earth standard).
- **Initial Angles (theta1, theta2):** Set to 1.2 and 2.1 radians, respectively.
- **Time Step (dt):** 0.0005 s, with `N_SUBSTEPS = 10` to update the screen at a more modest rate.
- **Simulation Duration:** 40.0 seconds (real time).
- **Frame Capture:** Optionally set `CAPTURE_FRAMES = True` to save each frame as a PNG.

**Feel free** to experiment with these parameters by editing their values in the script.

## References and Resources

- **Double Pendulum (Mathematical Background):**
  - T. L. Reber, *Chaotic Dynamics of the Double Pendulum*, American Journal of Physics, 2024.
  - R. Fitzpatrick, *Analysis of a Double Pendulum*, University of Texas, 2018.  
- **VPython Documentation:** [https://vpython.org](https://vpython.org)
- **Runge-Kutta Methods:** E. Hairer, S. P. Nørsett, G. Wanner, *Solving Ordinary Differential Equations I*, Springer, 2nd ed., 2009.

## Contact

For questions, feedback, or collaboration inquiries regarding this project:

**Param Somane**  
Email: [psomane@ucsd.edu](mailto:psomane@ucsd.edu)  
**CSE 291D Physics Simulation – Spring 2025**  
University of California, San Diego

## License

This project is licensed under the terms of the **MIT License**. See [LICENSE](LICENSE) for details.
