from math import pi


class EmugBase(object):
    def __init__(self):
        pass

    @property
    def myu(self):
        myu = 4*pi*1e-7
        return myu

    @property
    def eps(self):
        eps = 8.8541878128 * 1e-12
        return eps
