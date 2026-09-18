# FEM_RoomTempControlProject
Deliverable: A MATLAB finite element model of transient heat transfer in a room with controlled heating along the top wall.


Completion: This project was completed April 24th, 2026, for the University of Arizona course AME 431 - Numerical Methods in Fluid Mechanics and Heat Transfer.

Objective: The objective of the project was to develop a finite element method code to solve the two-dimensional heat conduction problem and use it to simulate temperature variation in a room, and asses the effectiveness of a deadbeat control method in maintaining a specified temperature at a set location inside the room.
[Project2-FEM.pdf](https://github.com/user-attachments/files/32399934/Project2-FEM.pdf)

Part 1: The first part of this project was to develop a baseline code to validate steady state results that were provided, so as to confirm said code. project_2_baseline_copy.m provides this code.

Part 2: The second part of this project was to then modify the baseline code in order to implement time dependent Dirichlet boundary conditions, and then verify the results compared to a given temperature history.

Part 3: The last part of this project was to then implement temperature control through a deadbeat approach through the top wall, at a set location within the room. Using this code, an analysis was performed on the effectiveness of the control method under varying conditions, utilizing the underlying FEM concept to explain the obtained results.
