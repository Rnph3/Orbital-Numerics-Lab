import math
import unittest
from orbitlab.physics import (State, acceleration, circular_initial, eccentric_initial,
    euler_step, verlet_step, rk4_step, integrate, diagnostics, energy,angular_momentum)

class OrbitTests(unittest.TestCase):
    def test_circular_energy(self): self.assertAlmostEqual(energy(circular_initial()),-.5)
    def test_circular_h(self): self.assertAlmostEqual(angular_momentum(circular_initial()),1.)
    def test_eccentric_initial_r(self):self.assertAlmostEqual(eccentric_initial(.4).x,.6)
    def test_eccentric_energy(self):self.assertAlmostEqual(energy(eccentric_initial(.4)),-.5)
    def test_invalid_eccentricity(self):
        for x in (-.1,1,1.5):
            with self.assertRaises(ValueError):eccentric_initial(x)
    def test_acceleration(self):self.assertAlmostEqual(acceleration(1,0)[0],-1.)
    def test_acceleration_origin(self):
        with self.assertRaises(ValueError):acceleration(0,0)
    def test_bad_settings(self):
        for kwargs in (dict(dt=0,steps=10),dict(dt=.1,steps=0),dict(dt=.1,steps=10,method='no')):
            with self.assertRaises(ValueError):integrate(circular_initial(),**kwargs)
    def test_length(self):self.assertEqual(len(integrate(circular_initial(),dt=.1,steps=7)),8)
    def test_short_step(self):
        s=circular_initial()
        for fn in (euler_step,verlet_step,rk4_step):self.assertLess(fn(s,1e-5).y,1e-4)
    def test_rk4_close_to_circle(self):
        s=integrate(circular_initial(),dt=2*math.pi/100,steps=100,method='rk4')
        self.assertLess(diagnostics(s)['endpoint_position_error'],1e-4)
    def test_verlet_energy(self):
        s=integrate(circular_initial(),dt=2*math.pi/100,steps=100,method='verlet')
        self.assertLess(diagnostics(s)['max_relative_energy_drift'],1e-4)
    def test_verlet_h(self):
        s=integrate(circular_initial(),dt=.05,steps=100,method='verlet')
        self.assertLess(diagnostics(s)['max_relative_angular_momentum_drift'],1e-12)
    def test_euler_worse_than_rk4(self):
        initial=circular_initial();dt=2*math.pi/100
        eu=diagnostics(integrate(initial,dt=dt,steps=100,method='euler'))
        rk=diagnostics(integrate(initial,dt=dt,steps=100,method='rk4'))
        self.assertGreater(eu['max_relative_energy_drift'],rk['max_relative_energy_drift']*1e3)
    def test_refinement_rk4(self):
        p=lambda n:diagnostics(integrate(circular_initial(),dt=2*math.pi/n,steps=n,method='rk4'))['endpoint_position_error']
        self.assertGreater(p(50),p(100)*10)
    def test_time_reversal_verlet(self):
        init=circular_initial();f=integrate(init,dt=.01,steps=100,method='verlet')[-1]
        # Velocity reversal uses positive dt: T(S) = (q, -v)
        reversed_state=State(f.x,f.y,-f.vx,-f.vy)
        r=integrate(reversed_state,dt=.01,steps=100,method='verlet')[-1]
        self.assertLess(math.hypot(r.x-init.x,r.y-init.y),1e-11)
    def test_angular_momentum_eccentric(self):
        s=integrate(eccentric_initial(.3),dt=.01,steps=100,method='verlet')
        self.assertLess(diagnostics(s)['max_relative_angular_momentum_drift'],1e-12)
    def test_invalid_diagnostics(self):
        with self.assertRaises(ValueError):diagnostics([])

if __name__=='__main__':unittest.main()
