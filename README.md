# Adaptive Control

Design and simulation of an **adaptive controller** that makes a two-degree-of-freedom system track a reference model, even though its parameters are unknown.
Course project for *Smart / Adaptive Control Systems*, Aristotle University of Thessaloniki (2024).

## What it does
- Defines the plant, a reference model and an adaptive control law with parameter-update laws
- Simulates the closed-loop system with `ode45`
- Plots how each state (θ₁, θ₂ and their velocities) tracks the reference model, and the tracking error over time

## Tools
MATLAB (`ode45`)

## How to run
Open MATLAB in this folder and run `EPSAE2024.m`.

## Report
Design, stability analysis and results: `EPSAE2024.pdf`.
