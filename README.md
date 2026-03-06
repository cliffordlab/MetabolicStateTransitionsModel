# MetabolicStateTransitionsModel
Code from paper: Modeling Metabolic State Transitions in Obesity Using a Time-Varying Lambda–Omega Framework


# Lambda–Omega Dynamics: Time-Varying Stability and Limit-Cycle Adaptation

This repository contains MATLAB scripts for simulating and visualizing **time-varying lambda–omega dynamical systems** under three different scenarios. The code generates animated phase-plane trajectories together with a time-series summary of peak amplitudes and cumulative peak-to-peak changes, and can optionally export each simulation as an `.mp4` movie.

The main goal of these scripts is to illustrate how changes in the radial growth term affect the system’s stability structure, nullclines, and limit-cycle behavior over time.

---

## Repository Overview

The repository includes three main simulation cases:

- `case_1.m` — **sign-flip transition** from stable oscillation to decay
- `case_2.m` — **growing limit cycle**
- `case_3.m` — **shrinking limit cycle**

Each script:
- solves a lambda–omega system using `ode45`
- animates the trajectory in the phase plane
- updates the vector field and nullclines over time
- tracks peaks of \(x(t)\)
- computes the cumulative absolute change between consecutive peaks
- highlights an adaptation window based on peak differences
- optionally saves the animation as a video

---

## Mathematical Background

The simulations are based on the planar lambda–omega system

\[
\dot{x} = \Lambda(r,t)x - \Omega(r,t)y
\]

\[
\dot{y} = \Omega(r,t)x + \Lambda(r,t)y
\]

where

\[
r = \sqrt{x^2 + y^2}
\]

Here:
- \(\Lambda(r,t)\) controls radial growth or decay
- \(\Omega(r,t)\) controls angular velocity

Depending on how \(\Lambda(r,t)\) is defined, the system can exhibit:
- a stable limit cycle
- a shrinking oscillation
- a growing oscillation
- a transition from oscillatory to non-oscillatory behavior

---

## Visualization Layout

Each simulation uses a two-panel figure:

### Left panel: Phase plane
Shows:
- trajectory in \((x,y)\)
- moving current state
- dynamic nullclines
- frozen nullcline snapshots at selected intervals
- vector field
- background shading indicating radial stability

### Right panel: Peak-based summary
Shows:
- peak values of \(x(t)\)
- cumulative absolute differences between consecutive peaks
- shaded adaptation region based on peak-change thresholding

---

## Case Descriptions

## Case 1 — Stability Flip / Oscillation Collapse

**File:** `case_1.m`

This case uses a lambda–omega model with a **time-varying sign change** in the radial term:

\[
\Lambda(r,t) = s(t)\lambda_{\text{base}}(r), \quad \lambda_{\text{base}}(r)=1-r^2
\]

where \(s(t)\) transitions smoothly from \(+1\) to \(-1\) using a cosine function.

### Interpretation
- Initially, the system supports a stable oscillation
- Over time, the sign flip reverses the radial stability
- The oscillation weakens and eventually collapses toward the origin

### What this case demonstrates
- a smooth transition in system stability
- deformation of nullclines during the transition
- reduction in oscillation amplitude over time
- detection of the adaptation period from peak differences

### Output
- animated phase portrait
- peak evolution over time
- cumulative peak-change curve
- optional movie export to:
  - `Movie/case_1_movie.mp4`

---

## Case 2 — Growing Limit Cycle

**File:** `case_2.m`

This case uses an order-2 lambda–omega system with:

\[
\Lambda(r,t)=\lambda(t)-br^2
\]

\[
\Omega(r,t)=\omega_0 + ar^2
\]

The target limit-cycle radius grows smoothly from a smaller initial radius to a larger final radius.

### Interpretation
- the oscillation begins with a smaller radius
- the target radius increases over time
- the trajectory expands outward until it approaches the larger limit cycle

### What this case demonstrates
- gradual outward growth of the oscillation
- time-varying nullclines associated with expanding amplitude
- increasing or adapting peak structure during the transient phase
- stabilization onto a larger oscillatory orbit

### Output
- animated phase portrait
- peak evolution over time
- cumulative peak-change curve
- optional movie export to:
  - `Movie/case_2_movie.mp4`

---

## Case 3 — Shrinking Limit Cycle

**File:** `case_3.m`

This case uses the same order-2 lambda–omega structure as Case 2, but the target limit-cycle radius now **shrinks smoothly** over time.

### Interpretation
- the oscillation starts on a larger initial cycle
- the target radius decreases gradually
- the trajectory contracts inward until it approaches a smaller limit cycle

### What this case demonstrates
- gradual inward contraction of oscillatory behavior
- time-varying nullcline motion during shrinking dynamics
- transient adaptation visible in peak amplitudes
- stabilization onto a smaller oscillatory orbit

### Output
- animated phase portrait
- peak evolution over time
- cumulative peak-change curve
- optional movie export to:
  - `Movie/case_3_movie.mp4`

---

## Requirements

This code was written in **MATLAB** and relies on standard MATLAB functions, including:

- `ode45`
- `odeset`
- `findpeaks`
- `tiledlayout`
- `nexttile`
- `quiver`
- `contour`
- `imagesc`
- `patch`
- `VideoWriter`

### MATLAB toolboxes
You may need:
- **MATLAB**
- **Signal Processing Toolbox** (for `findpeaks`)

---

## How to Run

1. Clone the repository:
   ```bash
   git clone <your-repository-url>
   cd <your-repository-folder>
