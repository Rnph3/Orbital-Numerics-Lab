"""A two-body test-particle orbit in nondimensional central gravity.

No N-body gravitational interaction, atmospheric drag, thrust, perturbations,
relativity, external ephemerides or mission flight dynamics are modeled.
Units: gravitational parameter mu=1; circular radius=1, speed=1, period=2*pi.
"""
from __future__ import annotations
from dataclasses import dataclass
from math import hypot, sqrt


@dataclass(frozen=True)
class State:
    x: float
    y: float
    vx: float
    vy: float


def acceleration(x: float, y: float, mu: float = 1.0) -> tuple[float, float]:
    r = hypot(x, y)
    if r <= 0 or mu <= 0:
        raise ValueError('Radius and gravitational parameter must be positive')
    factor = -mu / (r ** 3)
    return factor * x, factor * y


def derivative(s: State, mu: float = 1.0) -> State:
    ax, ay = acceleration(s.x, s.y, mu)
    return State(s.vx, s.vy, ax, ay)


def energy(s: State, mu: float = 1.0) -> float:
    return (s.vx ** 2 + s.vy ** 2) / 2 - mu / hypot(s.x, s.y)


def angular_momentum(s: State) -> float:
    return s.x * s.vy - s.y * s.vx


def circular_initial() -> State:
    return State(1., 0., 0., 1.)


def eccentric_initial(e: float) -> State:
    if not 0 <= e < 1:
        raise ValueError('Eccentricity must be in [0,1)')
    return State(1-e, 0., 0., sqrt((1+e)/(1-e)))


def add(s: State, d: State, factor: float) -> State:
    return State(s.x + factor*d.x, s.y + factor*d.y,
                 s.vx + factor*d.vx, s.vy + factor*d.vy)


def euler_step(s: State, dt: float, mu: float = 1.) -> State:
    return add(s, derivative(s, mu), dt)


def verlet_step(s: State, dt: float, mu: float = 1.) -> State:
    ax, ay = acceleration(s.x, s.y, mu)
    x2 = s.x + s.vx*dt + 0.5*ax*dt*dt
    y2 = s.y + s.vy*dt + 0.5*ay*dt*dt
    ax2, ay2 = acceleration(x2, y2, mu)
    return State(x2, y2, s.vx + 0.5*(ax+ax2)*dt, s.vy + 0.5*(ay+ay2)*dt)


def rk4_step(s: State, dt: float, mu: float = 1.) -> State:
    k1 = derivative(s,mu)
    k2 = derivative(add(s, k1, dt/2), mu)
    k3 = derivative(add(s, k2, dt/2), mu)
    k4 = derivative(add(s, k3, dt), mu)
    return State(*(getattr(s,f) + dt*(getattr(k1,f) + 2*getattr(k2,f) +
               2*getattr(k3,f) + getattr(k4,f))/6
                   for f in ('x','y','vx','vy')))


STEPPERS = {'euler':euler_step,'verlet':verlet_step,'rk4':rk4_step}


def integrate(initial:State, *, dt:float, steps:int, method:str='verlet', mu:float=1.) -> list[State]:
    if dt<=0 or steps<1 or method not in STEPPERS or mu<=0:
        raise ValueError('Invalid integrator settings')
    result=[initial]
    for _ in range(steps):
        result.append(STEPPERS[method](result[-1],dt,mu))
    return result


def diagnostics(states:list[State], mu:float=1.)->dict[str,float]:
    if not states:
        raise ValueError('No states provided')
    e0=energy(states[0],mu)
    h0=angular_momentum(states[0])
    if e0==0 or h0==0:
        raise ValueError('Undefined relative error reference')
    return {
       'max_relative_energy_drift':max(abs(energy(s,mu)-e0)/abs(e0) for s in states),
       'max_relative_angular_momentum_drift':max(abs(angular_momentum(s)-h0)/abs(h0) for s in states),
       'endpoint_position_error':hypot(states[-1].x-states[0].x,states[-1].y-states[0].y),
       'minimum_radius':min(hypot(s.x,s.y) for s in states),
       'maximum_radius':max(hypot(s.x,s.y) for s in states),
    }
