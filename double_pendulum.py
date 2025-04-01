#!/usr/bin/env python3
import math

import numpy as np
from vpython import (
    canvas, vector, color, sphere, cylinder,
    rate, distant_light, local_light, box
)

# ------------------------------------------------------------------------------
# 1) Helper function: Convert a hex color string to a normalized VPython vector
# ------------------------------------------------------------------------------
def hex_to_rgbnorm(hex_str):
    """Convert a hex color string (e.g. '#f26419') to a normalized RGB vector."""
    hex_str = hex_str.lstrip('#')  # remove leading '#' if present
    r = int(hex_str[0:2], 16) / 255.0
    g = int(hex_str[2:4], 16) / 255.0
    b = int(hex_str[4:6], 16) / 255.0
    return vector(r, g, b)

# ------------------------------------------------------------------------------
# 2) Define a new pastel color palette
# ------------------------------------------------------------------------------
# New palette: pastel pinks, lavenders, violets, and purples
BACKGROUND_PASTEL_HEX = '#FFF0F5'   # Lavender blush for a light background
PIVOT_COLOR_HEX       = '#EE82EE'   # Violet
ROD1_COLOR_HEX        = '#8F8788'   # Pastel pink
ROD2_COLOR_HEX        = '#8F8788'   # Pastel pink
MASS1_COLOR_HEX       = '#9370DB'   # Medium purple
MASS2_COLOR_HEX       = '#DA70D6'   # Orchid
MASS2_TRAIL_COLOR_HEX = '#DA70D6'   # Light orchid
FLOOR_COLOR_HEX       = '#D3D3D3'   # Light gray for the floor

# ------------------------------------------------------------------------------
# 3) Main double pendulum ODE + RK4 code
# ------------------------------------------------------------------------------
def double_pendulum_derivs(y, params):
    """
    Return time derivatives for a planar double pendulum.

    y = [theta1, omega1, theta2, omega2]
    params = (m1, m2, L1, L2, g)
    """
    theta1, omega1, theta2, omega2 = y
    m1, m2, L1, L2, g = params

    sin1 = np.sin(theta1)
    sin2 = np.sin(theta2)
    sin12 = np.sin(theta1 - theta2)
    cos12 = np.cos(theta1 - theta2)

    # Common denominator in the equations
    denom = 2*m1 + m2 - m2 * np.cos(2*theta1 - 2*theta2)

    # Angular acceleration for pendulum 1 (alpha1)
    alpha1 = (
        -g * (2*m1 + m2) * sin1
        - m2 * g * np.sin(theta1 - 2*theta2)
        - 2 * sin12 * m2 * (omega2**2 * L2 + omega1**2 * L1 * cos12)
    ) / (L1 * denom)

    # Angular acceleration for pendulum 2 (alpha2)
    alpha2 = (
        2 * sin12 *
        (
            omega1**2 * L1 * (m1 + m2)
            + g * (m1 + m2) * np.cos(theta1)
            + omega2**2 * L2 * m2 * cos12
        )
    ) / (L2 * denom)

    return np.array([omega1, alpha1, omega2, alpha2], dtype=float)

def rk4_step(y, dt, derivs_func, params):
    """Perform one 4th-order Runge-Kutta integration step."""
    k1 = derivs_func(y, params)
    k2 = derivs_func(y + 0.5*dt*k1, params)
    k3 = derivs_func(y + 0.5*dt*k2, params)
    k4 = derivs_func(y + dt*k3, params)
    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)

def run_simulation():
    # --- Physical / Simulation parameters ---
    m1, m2 = 2.0, 1.0
    L1, L2 = 1.5, 1.5
    g = 9.81
    params = (m1, m2, L1, L2, g)

    # Initial conditions
    theta1_0, omega1_0 = 1.2, 0.0
    theta2_0, omega2_0 = 2.1, 0.0
    y = np.array([theta1_0, omega1_0, theta2_0, omega2_0], dtype=float)

    dt = 0.0005
    sim_duration = 40.0

    # Number of small integration steps per *one* visual update
    N_SUBSTEPS = 10

    # Adjust 'steps' so total time is still 40s
    steps = int(sim_duration / dt / N_SUBSTEPS)

    # --- 3D scene setup (unchanged) ---
    scene = canvas(
        width=1280,
        height=720,
        center=vector(0, 1.0, 0),
        background=hex_to_rgbnorm(BACKGROUND_PASTEL_HEX),
        fov=0.0075
    )

    # Remove default lights; add custom lights
    scene.lights = []
    distant_light(direction=vector(1, 1, 1), color=color.white)
    local_light(pos=vector(-2, 3, 2), color=vector(0.7, 0.7, 0.7))
    local_light(pos=vector(2, -3, 3), color=vector(0.5, 0.5, 0.5))

    floor = box(
        pos=vector(0, -2.5, 0),
        size=vector(5, 0.05, 5),
        color=hex_to_rgbnorm(FLOOR_COLOR_HEX),
        opacity=0.2
    )

    pivot = vector(0, 1.5, 0)
    pivot_sphere = sphere(
        pos=pivot,
        radius=0.015,
        color=hex_to_rgbnorm(PIVOT_COLOR_HEX),
        shininess=0.6
    )

    rod1 = cylinder(
        pos=pivot,
        axis=vector(0, 0, 0),
        radius=0.03,
        color=hex_to_rgbnorm(ROD1_COLOR_HEX),
        shininess=0.6
    )
    rod2 = cylinder(
        pos=pivot,
        axis=vector(0, 0, 0),
        radius=0.03,
        color=hex_to_rgbnorm(ROD2_COLOR_HEX),
        shininess=0.6
    )

    mass1 = sphere(
        pos=vector(0, 0, 0),
        radius=0.1 * math.sqrt(2),
        color=hex_to_rgbnorm(MASS1_COLOR_HEX),
        shininess=0.6,
        make_trail=False
    )
    mass2 = sphere(
        pos=vector(0, 0, 0),
        radius=0.1,
        color=hex_to_rgbnorm(MASS2_COLOR_HEX),
        shininess=0.6,
        make_trail=True,
        trail_radius=0.01,
        trail_color=hex_to_rgbnorm(MASS2_TRAIL_COLOR_HEX),
        retain=5000
    )

    mass2.clear_trail()

    rod1.length = L1
    rod2.length = L2

    CAPTURE_FRAMES = False
    frame_count = 0

    # --- Main integration + rendering loop ---
    for i in range(steps):
        # 1) Do multiple small RK4 substeps to advance time
        for _ in range(N_SUBSTEPS):
            y = rk4_step(y, dt, double_pendulum_derivs, params)

        # 2) Update the display at a more modest rate
        rate(200)  # change this to 100 or 50 as you like

        theta1, omega1, theta2, omega2 = y

        # Compute new positions
        x1 = L1 * np.sin(theta1)
        y1 = -L1 * np.cos(theta1)
        x2 = x1 + L2 * np.sin(theta2)
        y2 = y1 - L2 * np.cos(theta2)

        # Shift so pivot is at (0, 1.5, 0)
        pos1 = vector(x1, y1, 0) + pivot
        pos2 = vector(x2, y2, 0) + pivot

        # Update rods + masses
        rod1.pos = pivot
        rod1.axis = pos1 - pivot
        rod2.pos = pos1
        rod2.axis = pos2 - pos1

        mass1.pos = pos1
        mass2.pos = pos2

        # Adjust the trail's opacity
        try:
            mass2.trail_object.opacity = 0.5
        except Exception:
            pass

        if CAPTURE_FRAMES:
            scene.capture(f"frame{i:04d}.png")
            frame_count += 1

    print("Simulation complete.")

if __name__ == "__main__":
    run_simulation()
