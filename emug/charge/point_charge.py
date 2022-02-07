from math import pi
import numpy as np

from emug import EmugBase


class PointCharge(EmugBase):
    def __init__(self, q):
        super(PointCharge, self).__init__()
        """
        Point charge Q positioned at (x,y,z)=(0,0,0)
        """
        self.q = q
    
    def ex(self, x, y, z):
        return self.q / (4*pi*self.eps) * x/(x**2+y**2+z**2)**(3/2)

    def ey(self, x, y, z):
        return self.q / (4*pi*self.eps) * y/(x**2+y**2+z**2)**(3/2)

    def ez(self, x, y, z):
        return self.q / (4*pi*self.eps) * z/(x**2+y**2+z**2)**(3/2)

    def e(self, x, y, z):
        return np.array([self.ex(x,y,z), self.ey(x,y,z), self.ez(x,y,z)])
        