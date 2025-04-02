#!/usr/bin/env python3
import math
import numpy as np
from vpython import (
    canvas, vector, color, sphere, cylinder, box,
    rate, distant_light, local_light
)

##############################################################################
# 1) Helper functions for color and small “indicator” functions
##############################################################################
def hex_to_rgbnorm(hex_str):
    """Convert a hex color string (e.g. '#f26419') to a normalized VPython vector."""
    hex_str = hex_str.lstrip('#')
    r = int(hex_str[0:2], 16) / 255.0
    g = int(hex_str[2:4], 16) / 255.0
    b = int(hex_str[4:6], 16) / 255.0
    return vector(r, g, b)

def sigma(j, k):
    """
    The sigma_{j,k} function from the reference:
      sigma(j,k) = 1 if j <= k, else 0.
    Note: j, k are assumed 0-based internally, so adjust carefully if needed.
    """
    return 1 if (j <= k) else 0

def phi(j, k):
    """
    The phi_{j,k} function from the reference:
      phi(j,k) = 1 if j != k, else 0.
    j, k are 0-based inside Python, but the math formula was for 1-based.
    The definition doesn't change for 0-based though: same logic.
    """
    return 0 if (j == k) else 1

##############################################################################
# 2) Build the mass matrix M and forcing vector f for the N-pendulum
#    using the big summation approach from e.g. Equation (86) in the text.
#
#    We'll assume:
#       angles = [theta_0, ..., theta_{N-1}]
#       omegas = [omega_0, ..., omega_{N-1}]
#       lengths = [L_0, ..., L_{N-1}]
#       masses  = [m_0, ..., m_{N-1}]
#       g = 9.81 (gravity)
#
#    The final system is:
#       M( angles, omegas ) * [alpha_0, ..., alpha_{N-1}]^T = - f( angles, omegas )
#
#    We will return (M, f). Then alpha = - np.linalg.solve(M, f).
##############################################################################

def build_M_and_f(angles, omegas, lengths, masses, g):
    """
    Construct the NxN matrix M and the N-vector f for the N-pendulum system.

    angles, omegas:  length N
    lengths, masses: length N
    """
    N = len(angles)
    Mmat = np.zeros((N, N), dtype=float)
    fvec = np.zeros(N, dtype=float)

    # For convenience of notation:
    th = angles
    w  = omegas
    L  = lengths
    m  = masses

    # We'll follow a direct summation approach reminiscent of eq. (86):
    #
    #   0 = sum_{k=1..N} [ g * L_j * sin(th_j)* m_k * sigma_{j,k}
    #                     + m_k * L_j^2 * ddot{th_j} * sigma_{j,k}
    #                     + (sum_{q >= k} m_q sigma_{j,q}) * L_j * L_k
    #                        * sin(th_j - th_k) * w_j^2
    #                     + ...
    #                     +  (some terms that multiply ddot{th_k})
    #                   ]
    #
    # We'll separate out everything that multiplies ddot{th_j} or ddot{th_k}
    # (which go into Mmat) from everything else (which goes to fvec).

    for j in range(N):
        # We'll accumulate eqn_j in the form:
        #   eqn_j = sum_of_stuff = 0
        # eventually we want eqn_j = Mmat[j,:]*alpha + fvec[j] = 0
        # => Mmat[j,:]*alpha = -fvec[j].
        eqn_j_constant = 0.0  # accumulates terms that do not multiply any ddot{th}

        #--------------------------------------------
        # (1) The term: g * L_j * sin(th_j) * sum_{k} [m_k sigma_{j,k}]
        #     does NOT multiply ddot{th} => goes to fvec
        #--------------------------------------------
        for k in range(N):
            s_jk = sigma(j, k)
            if s_jk == 1:
                eqn_j_constant += g * L[j] * math.sin(th[j]) * m[k]

        #--------------------------------------------
        # (2) The term: m_k * L_j^2 * ddot{th_j} * sigma_{j,k}
        #     *does* multiply ddot{th_j} => add to Mmat[j,j]
        #--------------------------------------------
        for k in range(N):
            s_jk = sigma(j, k)
            if s_jk == 1:
                Mmat[j, j] += m[k] * (L[j]**2)

        #--------------------------------------------
        # (3) The big velocity-coupling term:
        #    (sum_{q >= k} m_q sigma_{j,q}) * L_j * L_k * sin(th_j - th_k) * w_j * w_k
        #    -- Actually in eqn(86), we have two sets of velocity terms:
        #       sin(th_j - th_k) * w_j * w_k  and also the piece
        #       sin(th_k - th_j)[(w_j - w_k)*w_k], etc.
        #
        #    We'll break it down carefully:
        #
        #    3a) The term
        #        (sum_{q>=k} m_q sigma_{j,q}) * L_j*L_k * sin(th_j - th_k)* w_j*w_k
        #        does NOT multiply ddot{th}. => goes to fvec
        #
        #    3b) The sub-term that multiplies ddot{th_k} => will be put in Mmat[j,k].
        #
        # In the paper's final eqn(86), the “4th sum” has structure:
        #   + (sum_{q>=k} m_q sigma_{j,q}) * L_j L_k [
        #         sin(th_k - th_j)* (w_j - w_k)* w_k
        #       + phi_{j,k} cos(th_j - th_k)* ddot{th_k}   <-- the part that multiplies ddot{th_k}
        #   ]
        #
        # We'll do it in pieces.
        #--------------------------------------------
        for k in range(N):
            # We will define an effective mass factor: Mfactor = sum_{q >= k} m_q sigma_{j,q}
            # i.e. sum_{q=k..N-1 if q >= k} m_q * sigma_{j,q}.
            # But be careful with indexing: q >= k means q in [k..N-1].
            Mfactor = 0.0
            for q in range(k, N):
                Mfactor += m[q]*sigma(j, q)

            #----------------------------------------
            # 3a) Terms that do NOT multiply ddot{th_*}
            #     sin(th_j - th_k) * w_j * w_k
            #     Also: sin(th_k - th_j)*(...) from the paper
            #----------------------------------------
            # The paper’s eqn(86) has:
            #   + (Mfactor)* L_j * L_k * sin(th_j - th_k)* w_j*w_k
            #   + (Mfactor)* L_j * L_k * sin(th_k - th_j)[(w_j - w_k)* w_k]
            # In code, sin(th_k - th_j) = -sin(th_j - th_k).
            #
            # Let’s piece them out:

            # 3a-i) (Mfactor)* L_j L_k sin(th_j - th_k)* w_j*w_k
            eqn_j_constant += Mfactor * L[j]*L[k] * math.sin(th[j] - th[k]) * (w[j]*w[k])

            # 3a-ii) (Mfactor)* L_j L_k * sin(th_k - th_j)* ( (w_j - w_k)*w_k )
            #        = (Mfactor)* L_j L_k * [ -sin(th_j - th_k) ] * [ (w_j - w_k)*w_k ]
            eqn_j_constant += Mfactor * L[j]*L[k] * (
                - math.sin(th[j] - th[k]) * ((w[j] - w[k])*w[k])
            )

            #----------------------------------------
            # 3b) Terms that *do* multiply ddot{th_k}:
            #   (Mfactor)* L_j L_k * phi_{j,k} cos(th_j - th_k)* ddot{th_k}
            #
            # So we add to Mmat[j, k]:
            #     += (Mfactor)* L_j L_k * phi_{j,k} * cos(th_j - th_k)
            #----------------------------------------
            # Because ddot{th_k} is the unknown for mass k in row j's equation.
            Mmat[j, k] += Mfactor * L[j]*L[k] * phi(j, k) * math.cos(th[j] - th[k])

        #--------------------------------------------
        # Summation done.  eqn_j = 0 => means:
        #   sum_of_all_ddot_terms (which we put into Mmat) + eqn_j_constant = 0.
        # So fvec[j] = eqn_j_constant.  Then we solve Mmat * alpha = -fvec => alpha = -M^-1 f.
        #--------------------------------------------
        fvec[j] = eqn_j_constant

    return (Mmat, fvec)


##############################################################################
# 3) The ODE function: given y = [theta_0,omega_0, ..., theta_{N-1},omega_{N-1}],
#    return dydt = [omega_0, alpha_0, ..., omega_{N-1}, alpha_{N-1}].
##############################################################################

def n_pendulum_derivs(y, lengths, masses, g):
    """
    For an N-pendulum state vector y of length 2N, compute d/dt of y.
    We:
      1) parse angles & omegas
      2) build M, f
      3) solve M * alpha = -f for alpha
      4) return [omega_i, alpha_i]
    """
    N = len(masses)
    thetas = y[0:N]
    omegas = y[N:2*N]

    # Build the system:
    Mmat, fvec = build_M_and_f(thetas, omegas, lengths, masses, g)

    # Solve for alpha = - M^-1 f:
    # We do a linear solve each step.  If M is singular or ill-conditioned for certain
    # angles, that's a physical singularity (e.g. rods collinear).  For typical angles,
    # it's solvable.
    alpha = -np.linalg.solve(Mmat, fvec)

    # Construct the derivative array dydt:
    dydt = np.zeros_like(y)
    # first N are dtheta/dt = omega
    for i in range(N):
        dydt[i] = omegas[i]
    # next N are domega/dt = alpha
    for i in range(N):
        dydt[N+i] = alpha[i]

    return dydt


##############################################################################
# 4) 4th-order Runge–Kutta step
##############################################################################
def rk4_step(y, dt, derivs_func, *extra_args):
    """
    Perform one 4th-order Runge-Kutta integration step.
    derivs_func(y, *extra_args) -> derivative of y.
    """
    k1 = derivs_func(y, *extra_args)
    k2 = derivs_func(y + 0.5*dt*k1, *extra_args)
    k3 = derivs_func(y + 0.5*dt*k2, *extra_args)
    k4 = derivs_func(y + dt*k3, *extra_args)
    return y + (dt/6.0)*(k1 + 2*k2 + 2*k3 + k4)


##############################################################################
# 5) Main simulation routine
##############################################################################

def run_simulation():
    # -----------------------------
    # A) USER PARAMETERS
    # -----------------------------
    N = 8  # <--- Number of bobs (change this as desired)

    # Example: set all rods the same length, all masses the same,
    # or customize them as arrays of length N:
    # For demonstration, let's do random-ish or varied lengths and masses:
    lengths = np.array([1.5]*N, dtype=float)  # all rods ~1.5
    masses  = np.array([1.0]*N, dtype=float)  # all masses ~1.0
    g = 9.81

    # Initial angles (in radians) and angular velocities:
    # (Here we just pick a small spread of angles, all zero velocities.)
    # You MUST provide N angles and N angular velocities to fill the 2N state.
    init_thetas = np.array([0.5*(i+1) for i in range(N)], dtype=float)   # e.g. [0.3, 0.6, 0.9, 1.2,...]
    init_omegas = np.zeros(N, dtype=float)

    # Build initial state vector y0
    y0 = np.concatenate([init_thetas, init_omegas])

    dt = 0.001        # Timestep for RK4
    sim_duration = 30.0
    # Number of substeps per one visual update
    N_SUBSTEPS = 10
    steps = int(sim_duration / dt / N_SUBSTEPS)

    # -----------------------------
    # B) VPython 3D scene setup
    # -----------------------------
    # Feel free to customize colors:
    BACKGROUND_COLOR = hex_to_rgbnorm("#FFF0F5")  # lavender-ish
    ROD_COLOR        = hex_to_rgbnorm("#8F8788")
    FLOOR_COLOR      = hex_to_rgbnorm("#D3D3D3")
    # We'll pick a mild color gradient for the masses:
    MASS_COLORS_HEX  = ["#9370DB", "#DA70D6", "#EE82EE", "#F08080", "#FFA07A",
                        "#FFB6C1", "#FFA500", "#FA8072", "#BA55D3", "#7B68EE"]

    scene = canvas(width=1280, height=720, background=BACKGROUND_COLOR)
    scene.lights = []
    distant_light(direction=vector(1, 1, 1), color=color.white)
    local_light(pos=vector(-2, 3, 2), color=vector(0.7, 0.7, 0.7))
    local_light(pos=vector(2, -3, 3), color=vector(0.5, 0.5, 0.5))

    # A floor just for context:
    floor = box(pos=vector(0, -2.5, 0), size=vector(5, 0.05, 5),
                color=FLOOR_COLOR, opacity=0.2)

    # Let's place the top pivot at pivot_pos
    pivot_pos = vector(0, 1.5, 0)

    # Build rods & masses in a loop
    rods = []
    masses_spheres = []
    for i in range(N):
        rod = cylinder(pos=pivot_pos if i == 0 else vector(0,0,0),
                       axis=vector(0,0,0),
                       radius=0.02,
                       color=ROD_COLOR,
                       shininess=0.6)
        rods.append(rod)

        sp = sphere(pos=pivot_pos,  # dummy init
                    radius=0.06 + 0.02*i,  # slightly bigger each time, for fun
                    color=hex_to_rgbnorm(MASS_COLORS_HEX[i % len(MASS_COLORS_HEX)]),
                    shininess=0.6,
                    make_trail=(i == N-1),   # only the last mass leaves a trail, optional
                    trail_radius=0.01,
                    retain=2000)

        if i == N - 1:
            sp.clear_trail()

        masses_spheres.append(sp)

    # -----------------------------
    # C) Time Integration
    # -----------------------------
    y = np.copy(y0)

    for step_i in range(steps):
        # Take N_SUBSTEPS small RK4 steps
        for _ in range(N_SUBSTEPS):
            y = rk4_step(y, dt, n_pendulum_derivs, lengths, masses, g)

        # Update VPython visualization ~ (rate)
        rate(100)  # adjust as desired

        # Unpack the angles
        thetas = y[0:N]

        # Compute (x_i, y_i) for each bob i,
        #   x_i = sum_{k=0..i} L_k sin(theta_k)
        #   y_i = - sum_{k=0..i} L_k cos(theta_k)
        # Then shift by pivot_pos in the 3D scene
        x_coords = np.zeros(N)
        y_coords = np.zeros(N)
        for i in range(N):
            # x_i,y_i is the position of mass i in the plane (z=0)
            x_sum = 0.0
            y_sum = 0.0
            for k in range(i+1):
                x_sum += lengths[k]*math.sin(thetas[k])
                y_sum += lengths[k]*(-math.cos(thetas[k]))
            x_coords[i] = x_sum
            y_coords[i] = y_sum

        # Update rods & masses
        for i in range(N):
            # Mass i is at:
            pos_i = pivot_pos + vector(x_coords[i], y_coords[i], 0)

            masses_spheres[i].pos = pos_i

            # The rod i goes from:
            #   pivot_pos if i==0
            #   else the previous bob's position
            if i == 0:
                rods[i].pos = pivot_pos
            else:
                rods[i].pos = masses_spheres[i-1].pos

            rods[i].axis = masses_spheres[i].pos - rods[i].pos
            rods[i].length = rods[i].axis.mag

    print("Simulation complete.")

if __name__ == "__main__":
    run_simulation()
