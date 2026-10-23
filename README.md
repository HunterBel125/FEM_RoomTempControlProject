# Finite Element Room Temperature Control

**Deliverable:** A MATLAB finite element model of transient heat conduction in a two-dimensional room domain, with temperature control applied through the top-wall boundary.

**Completion:** This project was completed on April 24, 2026, for the University of Arizona course **AME 431 — Numerical Methods in Fluid Mechanics and Heat Transfer**.

**Objective:** Develop a finite element code to simulate temperature variation within a room domain and evaluate a deadbeat control method for maintaining a specified temperature at a selected interior location. The project progresses from steady-state verification to transient simulation and temperature control.

[View the project assignment](https://github.com/user-attachments/files/32399934/Project2-FEM.pdf)

## Part 1: Baseline Model and Steady-State Verification

[`project_2_baseline_copy.m`](project_2_baseline_copy.m) establishes the finite element model for the two-dimensional heat conduction problem.

The baseline code calculates the temperature distribution under prescribed boundary conditions. Its steady-state results are compared with the reference results provided in the assignment to check the implementation before introducing time-dependent behavior.

This stage establishes the foundation for the transient and controlled simulations.

## Part 2: Transient Heat Transfer

The second part extends the baseline model to simulate temperature changes over time and introduces time-dependent Dirichlet boundary conditions.

These boundary conditions prescribe temperature as a function of time, allowing the model to represent changing thermal conditions along the room boundary. The calculated temperature history is compared with the reference history supplied in the assignment to check the transient implementation.

This stage connects the steady-state model to the time-dependent response needed for temperature control.

## Part 3: Deadbeat Temperature Control

The final part introduces a deadbeat control approach that adjusts the top-wall temperature to regulate temperature at a selected location inside the room.

The controlled model is used to investigate how effectively the target temperature is maintained under varying conditions. The analysis connects the observed response to the underlying finite element model and the way heat propagates from the controlled boundary to the interior measurement location.

This stage brings together the numerical model, transient response, and control method to evaluate temperature regulation within the simulated domain.

## Verification Approach

The project uses two reference comparisons before evaluating the controller:

- **Steady-state verification:** Compare the baseline temperature distribution with the supplied steady-state results.
- **Transient verification:** Compare the calculated temperature history with the supplied time-dependent reference results.

These comparisons help assess whether the numerical implementation reproduces the expected behavior. Controller performance is then evaluated using the temperature response at the selected interior location.

The reference comparisons are numerical checks; they do not establish agreement with measurements from a physical room.

## Project Organization

The three parts follow a sequence of increasing capability:

1. Establish and check the baseline heat conduction model.
2. Extend the model to time-dependent boundary conditions.
3. Apply boundary temperature control and evaluate its effectiveness.

Each stage builds on the preceding model, allowing the thermal response to be checked before introducing the additional behavior required for control.

## Running the Baseline Model

Download the project files and set the project folder as MATLAB’s Current Folder. Run the baseline script with:

```matlab
project_2_baseline_copy
```

Review the generated results against the assignment’s reference case before proceeding to the transient and controlled versions.

## Scope and Interpretation

This project models heat conduction in a two-dimensional domain representing a room. The room interpretation provides a setting for studying numerical heat transfer and control; it should not be treated as a complete building thermal model.

The controlled input is the prescribed top-wall temperature. This represents an idealized thermal boundary condition, rather than a detailed model of a physical heater or HVAC system.

The main focus is the implementation and verification of the finite element method, followed by an assessment of how the modeled thermal response influences temperature-control performance.

## Skills Applied

MATLAB programming, finite element modeling, transient heat conduction, time-dependent boundary conditions, numerical verification, temperature control, and interpretation of simulation results.
