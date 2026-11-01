# A Correct Gravity Equation Can Still Produce the Wrong Orbit

*Comparing Euler, Velocity Verlet, and RK4 by energy drift, endpoint error, and grid refinement.*

A numerical simulation can obey the correct physical equation and still report nonphysical motion. The integration method and step size are part of the scientific claim, not incidental implementation details.

I built **Orbital Numerics Lab** to compare three time integrators using the same simple gravitational problem.

## A controlled Kepler experiment

The model uses a test particle in a fixed inverse-square central gravitational field. In normalized units, the gravitational parameter is μ = 1, the initial circular radius is 1, the speed is 1, and one orbital period is **2π**.

The acceleration is:

**a(r) = −μr / |r|³**

Two useful physical quantities are the specific energy, E = |v|²/2 − μ/|r|, and angular momentum, h = x·vy − y·vx. Both should remain constant in the underlying ideal two-body model.

## Three integrators, one time interval

- **Explicit Euler** advances with the current derivative. It is simple and only first-order accurate globally.
- **Velocity Verlet** uses acceleration at the start and end of a step, and has useful symplectic and time-reversible properties for this conservative system.
- **Classic RK4** evaluates four slopes per step and is fourth-order accurate globally under standard smoothness assumptions.

I used **160 equal steps over exactly one period** for every method. The program then measured maximum relative energy drift and the distance between the final position and the starting point. It also evaluated refinement from 40 to 320 steps per period.

## Results

| Method | Maximum relative energy drift | Endpoint position error |
|---|---:|---:|
| Explicit Euler | 2.6485 × 10⁻¹ | 1.6745 |
| Velocity Verlet | 5.9317 × 10⁻⁷ | 0.0032279 |
| RK4 | 1.6305 × 10⁻⁸ | 4.1917 × 10⁻⁷ |

[Insert `artifacts/orbit_comparison.png` and `artifacts/convergence.png` here.]

Under these **short-horizon, normalized conditions**, RK4 produced the smallest endpoint and energy errors. That does not mean Verlet is inferior in every situation: symplectic methods are especially important when long-time qualitative behavior matters, and a one-orbit test cannot establish a long-term advantage.

The Euler simulation, meanwhile, visibly drifts away from the intended orbit despite starting from the same correct equations.

## Why the engineering choices matter

The implementation separates physics (`acceleration`, `energy`, `angular_momentum`) from integration methods (`euler_step`, `verlet_step`, `rk4_step`) and from the experiment runner. That makes the numerical method independently testable. I also included a velocity-reversal test for Verlet and grid-refinement curves to avoid relying on a single resolution.

## Limits and next work

This is a **dimensionless two-body educational experiment**. No mission ephemeris, three-body perturbations, drag, J2, relativistic corrections, propulsion, state-estimation noise, or spacecraft control have been tested. Real mission-grade accuracy is far outside the result demonstrated here.

I would next compare thousands of orbits, multiple eccentricities, adaptive solvers, and an N-body extension. A useful frontier step would be connecting verified numerical dynamics with estimation and control algorithms.

**Code, 18 unit tests, charts, and methods:** [INSERT VERIFIED PUBLIC GITHUB REPOSITORY URL]

References: [Rutgers numerical errors notes](https://www.physics.rutgers.edu/~haule/ComputationalPhysics/Numerical_Error.html) · [Illinois numerical integrators module](https://heath.cs.illinois.edu/iem/ode/symplectic/).
