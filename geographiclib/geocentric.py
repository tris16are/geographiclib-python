"""geocentric.py: transcription of GeographicLib::Geocentric class."""
# geocentric.py
#
# This is a rather literal translation of the GeographicLib::Geocentric class
# to python.  See the documentation for the C++ class for more information at
#
#    https://geographiclib.sourceforge.io/C++/doc/annotated.html
#
# Copyright (c) Charles Karney (2008-2022) <karney@alum.mit.edu> and
# licensed under the MIT/X11 License.  For more information, see
# https://geographiclib.sourceforge.io/
######################################################################

from typing import TypeVar
import math
import sys

from geographiclib.geomath import Math
from geographiclib.constants import Constants


T = TypeVar('T', bound='Geocentric')

class Geocentric:
  """
  Convert between geodetic and geocentric coordinates.

  Convert between geodetic coordinates latitude = lat, longitude = lon,
  height = h (measured vertically from the surface of the ellipsoid) to
  geocentric coordinates (X, Y, Z).  The origin of geocentric coordinates
  is at the center of the earth.  The Z axis goes through the north pole,
  lat = 90°.  The X axis goes through lat = 0, lon = 0.  Geocentric
  coordinates are also known as earth centered, earth fixed (ECEF)
  coordinates.

  The conversion from geographic to geocentric coordinates is
  straightforward.  For the reverse transformation we use H. Vermeille,
  "Direct transformation from geocentric coordinates to geodetic
  coordinates", J. Geodesy 76, 451-454 (2002).

  Several changes have been made to ensure that the method returns accurate
  results for all finite inputs (even if h is infinite).  See C. F. F. Karney,
  "Geodesics on an ellipsoid of revolution", arxiv:1102.1215v1 (2011),
  Appendix B.
  """

  def __init__(self, a, f):
    """
    Construct a Geocentric object.

    :param a: the equatorial radius of the ellipsoid (meters)
    :param f: the flattening of the ellipsoid

    An exception is thrown if either of the axes of the ellipsoid is
    not positive.
    """
    self.a = float(a)
    """The equatorial radius in meters (readonly)"""
    self.f = float(f)
    """The flattening (readonly)"""
    self._e2 = self.f * (2 - self.f)
    self._e2m = Math.sq(1 - self.f)  # 1 - _e2
    self._e2a = abs(self._e2)
    self._e4a = Math.sq(self._e2)
    self._maxrad = 2 * self.a / sys.float_info.epsilon

    if not (math.isfinite(self.a) and self.a > 0):
      raise ValueError("Equatorial radius is not positive")
    if not (math.isfinite(self.f) and self.f < 1):
      raise ValueError("Polar semi-axis is not positive")

  @staticmethod
  def _Rotation(sphi, cphi, slam, clam):
    """
    Compute the rotation matrix from local ENU to geocentric XYZ coordinates.

    :param sphi: sine of latitude
    :param cphi: cosine of latitude
    :param slam: sine of longitude
    :param clam: cosine of longitude
    :return: 9-element list representing 3x3 rotation matrix in row-major order

    The rotation matrix M transforms a vector v in local ENU (East, North, Up)
    coordinates to geocentric XYZ coordinates: v_XYZ = M · v_ENU
    """
    # This rotation matrix is given by the following quaternion operations
    # qrot(lam, [0,0,1]) * qrot(phi, [0,-1,0]) * [1,1,1,1]/2
    # or
    # qrot(pi/2 + lam, [0,0,1]) * qrot(-pi/2 + phi , [-1,0,0])
    return [
        # Local X axis (east) in geocentric coords
        -slam,        clam,         0,
        # Local Y axis (north) in geocentric coords
        -clam * sphi, -slam * sphi, cphi,
        # Local Z axis (up) in geocentric coords
        clam * cphi,  slam * cphi,  sphi
    ]

  @staticmethod
  def _Rotate(M, x, y, z):
    """
    Apply rotation matrix M to transform from local to geocentric coordinates.

    :param M: 9-element rotation matrix in row-major order
    :param x: local x coordinate (east)
    :param y: local y coordinate (north)
    :param z: local z coordinate (up)
    :return: tuple (X, Y, Z) in geocentric coordinates

    Performs [X, Y, Z]^T = M · [x, y, z]^T
    """
    X = M[0] * x + M[1] * y + M[2] * z
    Y = M[3] * x + M[4] * y + M[5] * z
    Z = M[6] * x + M[7] * y + M[8] * z
    return X, Y, Z

  @staticmethod
  def _Unrotate(M, X, Y, Z):
    """
    Apply inverse rotation matrix M^T to transform from geocentric to local.

    :param M: 9-element rotation matrix in row-major order
    :param X: geocentric X coordinate
    :param Y: geocentric Y coordinate
    :param Z: geocentric Z coordinate
    :return: tuple (x, y, z) in local ENU coordinates

    Performs [x, y, z]^T = M^T · [X, Y, Z]^T
    """
    x = M[0] * X + M[3] * Y + M[6] * Z
    y = M[1] * X + M[4] * Y + M[7] * Z
    z = M[2] * X + M[5] * Y + M[8] * Z
    return x, y, z

  def Forward(self, lat, lon, h, M=False):
    """
    Convert from geodetic to geocentric coordinates.

    :param lat: latitude of point (degrees)
    :param lon: longitude of point (degrees)
    :param h: height of point above the ellipsoid (meters)
    :param M: if True, include the rotation matrix in the result
    :return: a dictionary with the following key/value pairs:

      * X: geocentric coordinate (meters)
      * Y: geocentric coordinate (meters)
      * Z: geocentric coordinate (meters)
      * M: (only if M=True) rotation matrix as 9-element list (row-major order)

    lat should be in the range [-90°, 90°].

    If M is True, the rotation matrix M transforms a unit vector v in local
    ENU (East, North, Up) coordinates to geocentric XYZ coordinates via
    v_XYZ = M · v_ENU.
    """
    sphi, cphi = Math.sincosd(Math.LatFix(lat))
    slam, clam = Math.sincosd(lon)
    n = self.a / math.sqrt(1 - self._e2 * Math.sq(sphi))
    Z = (self._e2m * n + h) * sphi
    X = (n + h) * cphi
    Y = X * slam
    X *= clam
    result = {'X': X, 'Y': Y, 'Z': Z}
    if M:
      result['M'] = self._Rotation(sphi, cphi, slam, clam)
    return result

  def Reverse(self, X, Y, Z, M=False):
    """
    Convert from geocentric to geodetic coordinates.

    :param X: geocentric coordinate (meters)
    :param Y: geocentric coordinate (meters)
    :param Z: geocentric coordinate (meters)
    :param M: if True, include the rotation matrix in the result
    :return: a dictionary with the following key/value pairs:

      * lat: latitude of point (degrees)
      * lon: longitude of point (degrees)
      * h: height of point above the ellipsoid (meters)
      * M: (only if M=True) rotation matrix as 9-element list (row-major order)

    In general, there are multiple solutions and the result which minimizes
    the absolute value of h is returned, i.e., (lat, lon) corresponds to the
    closest point on the ellipsoid.  The value of lon returned is in the
    range [-180°, 180°].

    If M is True, the rotation matrix M transforms a unit vector v in local
    ENU (East, North, Up) coordinates to geocentric XYZ coordinates via
    v_XYZ = M · v_ENU. The inverse transformation is v_ENU = M^T · v_XYZ.
    """
    R = math.hypot(X, Y)
    slam = Y / R if R != 0 else 0
    clam = X / R if R != 0 else 1
    h = math.hypot(R, Z)  # Distance to center of earth

    if h > self._maxrad:
      # We are really far away (> 12 million light years); treat the earth
      # as a point and h, above, is an acceptable approximation to the
      # height.  This avoids overflow, e.g., in the computation of disc
      # below.  It's possible that h has overflowed to inf; but that's OK.
      #
      # Treat the case X, Y finite, but R overflows to +inf by scaling by 2.
      R = math.hypot(X / 2, Y / 2)
      slam = (Y / 2) / R if R != 0 else 0
      clam = (X / 2) / R if R != 0 else 1
      H = math.hypot(Z / 2, R)
      sphi = (Z / 2) / H
      cphi = R / H
    elif self._e4a == 0:
      # Treat the spherical case.  Dealing with underflow in the general
      # case with _e2 = 0 is difficult.  Origin maps to N pole same as
      # with ellipsoid.
      H = math.hypot(1 if h == 0 else Z, R)
      sphi = (1 if h == 0 else Z) / H
      cphi = R / H
      h -= self.a
    else:
      # Treat prolate spheroids by swapping R and Z here and by switching
      # the arguments to phi = atan2(...) at the end.
      p = Math.sq(R / self.a)
      q = self._e2m * Math.sq(Z / self.a)
      r = (p + q - self._e4a) / 6
      if self.f < 0:
        p, q = q, p
      if not (self._e4a * q == 0 and r <= 0):
        # Avoid possible division by zero when r = 0 by multiplying
        # equations for s and t by r^3 and r, resp.
        S = self._e4a * p * q / 4  # S = r^3 * s
        r2 = Math.sq(r)
        r3 = r * r2
        disc = S * (2 * r3 + S)
        u = r
        if disc >= 0:
          T3 = S + r3
          # Pick the sign on the sqrt to maximize abs(T3).  This minimizes
          # loss of precision due to cancellation.  The result is unchanged
          # because of the way T is used in the definition of u.
          T3 += -math.sqrt(disc) if T3 < 0 else math.sqrt(disc)  # T3 = (r * t)^3
          # N.B. cbrt always returns the real root.  cbrt(-8) = -2.
          T = Math.cbrt(T3)  # T = r * t
          # T can be zero; but then r2 / T -> 0.
          u += T + (r2 / T if T != 0 else 0)
        else:
          # T is complex, but the way u is defined the result is real.
          ang = math.atan2(math.sqrt(-disc), -(S + r3))
          # There are three possible cube roots.  We choose the root which
          # avoids cancellation.  Note that disc < 0 implies that r < 0.
          u += 2 * r * math.cos(ang / 3)

        v = math.sqrt(Math.sq(u) + self._e4a * q)  # guaranteed positive
        # Avoid loss of accuracy when u < 0.  Underflow doesn't occur in
        # e4 * q / (v - u) because u ~ e^4 when q is small and u < 0.
        uv = self._e4a * q / (v - u) if u < 0 else u + v  # u+v, guaranteed positive
        # Need to guard against w going negative due to roundoff in uv - q.
        w = max(0.0, self._e2a * (uv - q) / (2 * v))
        # Rearrange expression for k to avoid loss of accuracy due to
        # subtraction.  Division by 0 not possible because uv > 0, w >= 0.
        k = uv / (math.sqrt(uv + Math.sq(w)) + w)
        k1 = k if self.f >= 0 else k - self._e2
        k2 = k + self._e2 if self.f >= 0 else k
        d = k1 * R / k2
        H = math.hypot(Z / k1, R / k2)
        sphi = (Z / k1) / H
        cphi = (R / k2) / H
        h = (1 - self._e2m / k1) * math.hypot(d, Z)
      else:  # e4 * q == 0 && r <= 0
        # This leads to k = 0 (oblate, equatorial plane) and k + e^2 = 0
        # (prolate, rotation axis) and the generation of 0/0 in the general
        # formulas for phi and h.  using the general formula and target_lon by 0
        # in formula for h.  So handle this case by taking the limits:
        # f > 0: z -> 0, k      ->   e2 * sqrt(q)/sqrt(e4 - p)
        # f < 0: R -> 0, k + e2 -> - e2 * sqrt(q)/sqrt(e4 - p)
        zz = math.sqrt((self._e4a - p if self.f >= 0 else p) / self._e2m)
        xx = math.sqrt(self._e4a - p if self.f < 0 else p)
        H = math.hypot(zz, xx)
        sphi = zz / H
        cphi = xx / H
        if Z < 0:
          sphi = -sphi  # for tiny negative Z (not for prolate)
        h = -self.a * (self._e2m if self.f >= 0 else 1) * H / self._e2a

    lat = Math.atan2d(sphi, cphi)
    lon = Math.atan2d(slam, clam)
    result = {'lat': lat, 'lon': lon, 'h': h}
    if M:
      result['M'] = self._Rotation(sphi, cphi, slam, clam)
    return result

  def EquatorialRadius(self):
    """Return the equatorial radius of the ellipsoid (meters)."""
    return self.a

  def Flattening(self):
    """Return the flattening of the ellipsoid."""
    return self.f

  @classmethod
  def WGS84(cls) -> T:
      """Instantiation for the WGS84 ellipsoid"""
      return cls(a=Constants.WGS84_a, f=Constants.WGS84_f)
