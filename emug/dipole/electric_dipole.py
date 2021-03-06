from math import pi
import numpy as np

from emug import EmugBase


class ElectricDipole(EmugBase):
    def __init__(self, moment, freq, dipole_direction=2):
        super(ElectricDipole, self).__init__()
        """
        原点に位置する微小電流素片
        """
        self.freq = freq
        self.moment = moment
        self.omega = 2*pi*self.freq
        self.k = self.omega*np.sqrt(self.eps*self.myu)
        self.dipole_direction = dipole_direction

    def _axis_transpose(self, x, y, z):
        if (self.dipole_direction == 'x') or (self.dipole_direction == 0):
            self.dipole_direction = 0
            x, z = z, -x
        elif (self.dipole_direction == 'y') or (self.dipole_direction == 1):
            self.dipole_direction = 1
            y, z = z, -y
        elif (self.dipole_direction == 'z') or (self.dipole_direction == 2):
            self.dipole_direction = 2
        else:
            raise AttributeError
        return x, y, z

    def _er(self, r, t, theta):
        jwt = (self.omega*t)*1j
        jkr = (self.k*r)*1j
        return np.sqrt(self.myu / self.eps) * self.moment * np.exp(jwt - jkr) / (2*pi*r*r) * (1 + 1/jkr) * np.cos(theta)

    def _etheta(self, r, t, theta):
        jwt = (self.omega*t)*1j
        jkr = (self.k*r)*1j
        kr = self.k*r
        imaginary_component = np.sqrt(self.myu/self.eps) * self.moment * np.exp(jwt - jkr) * self.k/(4*pi*r) * (1 + 1/jkr - 1/(kr*kr)) * np.sin(theta)
        return (imaginary_component)*1j

    def _ephai(self, r, t, theta):
        return np.zeros_like(r)

    def _hr(self, r, t, theta):
        return np.zeros_like(r)

    def _htheta(self, r, t, theta):
        return np.zeros_like(r)

    def _hphai(self, r, t, theta):
        jwt = (self.omega*t)*1j
        jkr = (self.k*r)*1j
        imaginary_component = self.moment * np.exp(jwt - jkr) * self.k/(4*pi*r) * (1 + 1/jkr) * np.sin(theta)
        return (imaginary_component)*1j

    def _ex(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        phai = np.arctan2(y, x)
        er_x = self._er(r=r, t=t, theta=theta) * np.sin(theta) * np.cos(phai)
        etheta_x = self._etheta(r=r, t=t, theta=theta) * np.cos(theta) * np.cos(phai)
        ephai_x = -1 * self._ephai(r=r, t=t, theta=theta) * np.sin(phai)
        return er_x + etheta_x + ephai_x

    def _ey(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        phai = np.arctan2(y, x)
        er_y = self._er(r=r, t=t, theta=theta) * np.sin(theta) * np.sin(phai)
        etheta_y = self._etheta(r=r, t=t, theta=theta) * np.cos(theta) * np.sin(phai)
        ephai_y = self._ephai(r=r, t=t, theta=theta) * np.cos(theta)
        return er_y + etheta_y + ephai_y

    def _ez(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        er_z = self._er(r=r, t=t, theta=theta) * np.cos(theta)
        etheta_z = -1 * self._etheta(r=r, t=t, theta=theta) * np.sin(theta)
        return er_z + etheta_z

    def _hx(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        phai = np.arctan2(y, x)
        hphai_x = -1 * self._hphai(r=r, t=t, theta=theta) * np.sin(phai)
        return hphai_x

    def _hy(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        phai = np.arctan2(y, x)
        hphai_y = self._hphai(r=r, t=t, theta=theta) * np.cos(phai)
        return hphai_y

    def _hz(self, x, y, z, t):
        r = np.sqrt(x*x + y*y + z*z)
        theta = np.arctan2(np.sqrt(x*x+y*y), z)
        hr_z = self._hr(r=r, t=t, theta=theta) * np.cos(theta)
        htheta_z = -1 * self._htheta(r=r, t=t, theta=theta) * np.sin(theta)
        return hr_z + htheta_z

    def e(self, x, y, z, t):
        x, y, z = self._axis_transpose(x, y, z)
        e = np.array([self._ex(x, y, z, t), self._ey(x, y, z, t), self._ez(x, y, z, t)])
        e[self.dipole_direction], e[2] = e[2], -e[self.dipole_direction]
        return e

    def h(self, x, y, z, t):
        x, y, z = self._axis_transpose(x, y, z)
        h = np.array([self._hx(x, y, z, t), self._hy(x, y, z, t), self._hz(x, y, z, t)])
        h[self.dipole_direction], h[2] = h[2], -h[self.dipole_direction]
        return h

    def em(self, x, y, z, t):
        return np.concatenate([self.e(x,y,z,t), self.h(x,y,z,t)], axis=0)
