from __future__ import annotations
import numpy as np
from typing import TypeVar, Type


class Interval:
    def __init__(self, values: np.ndarray):
        self.values = values

    def min(self):
        return self.values.min()

    def max(self):
        return self.values.max()

    def as_index_of(self, X):
        return X[self.values]


T = TypeVar("T", bound="Points")


class Points:
    def __init__(self, xy: np.ndarray):
        assert xy.ndim == 2, "xy should be a 2-dimensional array"
        self.xy = xy

    def __getitem__(self, key):
        return Points(self.xy[key])

    def __setitem__(self, key, value):
        self.xy[key] = value

    @classmethod
    def from_xy(cls: Type[T], x: np.ndarray, y: np.ndarray) -> T:
        return cls(xy=np.column_stack((x, y)))

    @property
    def x(self):
        return self.xy[:, 0]

    @property
    def y(self):
        return self.xy[:, 1]

    def __iter__(self):
        yield self.x
        yield self.y


class Curve(Points):
    def __init__(self, xy: np.ndarray):
        super().__init__(xy)

    def axis_values(self, axis):
        if axis not in ["x", "y"]:
            raise ValueError(f"Axis must be in: {axis}")
        return getattr(self, axis)

    def critical_points(self):
        first_derivative = self.dydx(order=1)
        return first_derivative.crosses_zero_at()

    def maxima(self):
        dens_crits = self.critical_points()
        second_derivative = self.dydx(order=2)
        return np.array(
            [xx for xx in dens_crits if second_derivative.y[xx] < 0], dtype=np.int64
        )

    def minima(self):
        dens_crits = self.critical_points()
        second_derivative = self.dydx(order=2)
        return np.array(
            [xx for xx in dens_crits if second_derivative.y[xx] > 0], dtype=np.int64
        )

    def dydx(self, order=1) -> Curve:
        if order == 0:
            return self
        else:
            dydx = np.gradient(self.y, self.x)

            return Curve.from_xy(self.x, dydx).dydx(order=(order - 1))

    def crosses_threshold_at(self, threshold, V_ignore_sign=True, axis="y"):
        V = self.axis_values(axis=axis)
        if V_ignore_sign:
            V = abs(V)
        return np.where(np.diff(np.sign(V - threshold)))[0]

    def crosses_zero_at(self, axis="y"):
        if axis not in ["x", "y"]:
            raise ValueError(f"Axis must be in: {axis}")
        return self.crosses_threshold_at(threshold=0.0, V_ignore_sign=False, axis=axis)

    def threshold_intervals(self, threshold, V_ignore_sign=True, axis="y"):
        V = self.axis_values(axis=axis)
        if threshold == 0.0:
            return [
                Interval(s)
                for s in np.split(range(0, len(V)), self.crosses_zero_at(axis=axis))
            ]
        else:
            return [
                Interval(s)
                for s in np.split(
                    range(0, len(V)),
                    self.crosses_threshold_at(
                        threshold, V_ignore_sign=V_ignore_sign, axis=axis
                    ),
                )
            ]

