# import modules
import numpy as np
from numpy.linalg import norm

from galaxy.timecourse import TimeCourse

G = 4.498768e-6 # kpc^3 / M_sun / Gyr

class M33AnalyticOrbit():
    """
    """

    def __init__(self, outfile='./m33_orbits.txt'):

        self.outfile = outfile
        self.tc = TimeCourse()

        # CoM position and velocity of M33 relative to M31
        M31_p, M31_v = self.tc.get_one_com('M31', 0)
        M33_p, M33_v = self.tc.get_one_com('M33', 0)
        self.pos = M33_p - M31_p
        self.vel = M33_v - M31_v

        # set component scale lengths and masses for M31
        self.rdisk = 5 # kpc
        self.mdisk = 12.014e10 # M_sun

        self.rbulge = 1 # kpc
        self.mbulge = 1.905e10 # M_sun

        self.rhalo = 60.1 # kpc, Hernquist a from Hwk 5
        self.mhalo = 192.1e10 # M_sun

    def hernquist_accel(self, M, r_a, pos):
        """
        Calculates acceleration induced by a a Hernquist profile.
        Relevant for halo and bulge.

        Args:
            M (float):
                mass (implicit M_sun)
            r_a (float):
                Herquist scale factor (implicit kpc)
            pos (3-vector):
                (x,y,z) position as np.array (implicit kpc)
        
        Returns:
            gravitational acceleration (3-vector)
        """

        r = norm(pos)
        factor = -G * M / (r * (r + r_a)**2)
        return factor * pos

    def miyamoto_nagai_accel(self, M, r_d, pos):
        """
        """

        x,y,z = pos
        R = np.sqrt(x**2 + y**2)
        z_d = self.rdisk / 5.
        B = r_d + np.sqrt(z**2 + z_d**2)
        factor_x = -G * M / (R**2 + B**2)**1.5
        factor_y = factor_x
        factor_z = -G * M * B / (R**2 + B**2)**1.5 / np.sqrt(z**2 + z_d**2)

        factor = np.array([factor_x, factor_y, factor_z])

        return factor * pos

    def m31_accel(self, pos):
        """
        """

        halo_acc = self.hernquist_accel(self.mhalo, self.rhalo, pos)
        bulge_acc = self.hernquist_accel(self.mbulge, self.rbulge, pos)
        disk_acc = self.miyamoto_nagai_accel(self.mdisk, self.rdisk, pos)

        return halo_acc + bulge_acc + disk_acc

    def leapfrog(self, dt, pos, vel):
        """
        """

        pos_mid = pos + vel * dt/2
        acc_mid = self.m31_accel(pos_mid)
        vel_next = vel + acc_mid * dt
        pos_next = pos + (vel + vel_next) * dt/2

        return pos_next, vel_next

    def orbit_integrator(self, dt, t_0=0, t_max=10):
        """
        """

        # starting parameters
        pos = self.pos
        vel = self.vel
        t = t_0

        # slight uncertainty about how many rows we end up with
        # so store results in a list for now
        orbits = [ (t, *tuple(pos), *tuple(vel)), ]

        # loop over the timesteps
        while t <= t_max:
            pos, vel = self.leapfrog(dt, pos, vel)
            orbits.append( (t, *tuple(pos), *tuple(vel) ) )
            t += dt

        orbits_array = np.array(orbits)
        np.savetxt(self.outfile, orbits_array, fmt = "%11.3f"*7, comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
