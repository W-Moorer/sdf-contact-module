import numpy as np


class Function:
    def value(self, t):
        raise NotImplementedError

    def first_derivative(self, t):
        raise NotImplementedError

    def second_derivative(self, t):
        raise NotImplementedError


class ConstantFunction(Function):
    def __init__(self, c):
        self.c = float(c)

    def value(self, t):
        return self.c

    def first_derivative(self, t):
        return 0.0

    def second_derivative(self, t):
        return 0.0


class LinearFunction(Function):
    def __init__(self, a, b):
        self.a = float(a)
        self.b = float(b)

    def value(self, t):
        return self.a + self.b * t

    def first_derivative(self, t):
        return self.b

    def second_derivative(self, t):
        return 0.0


class SinFunction(Function):
    def __init__(self, amplitude, omega, phase=0.0):
        self.A = float(amplitude)
        self.omega = float(omega)
        self.phase = float(phase)

    def value(self, t):
        return self.A * np.sin(self.omega * t + self.phase)

    def first_derivative(self, t):
        return self.A * self.omega * np.cos(self.omega * t + self.phase)

    def second_derivative(self, t):
        return -self.A * self.omega**2 * np.sin(self.omega * t + self.phase)


class CosFunction(Function):
    def __init__(self, amplitude, omega, phase=0.0):
        self.A = float(amplitude)
        self.omega = float(omega)
        self.phase = float(phase)

    def value(self, t):
        return self.A * np.cos(self.omega * t + self.phase)

    def first_derivative(self, t):
        return -self.A * self.omega * np.sin(self.omega * t + self.phase)

    def second_derivative(self, t):
        return -self.A * self.omega**2 * np.cos(self.omega * t + self.phase)
