"""Geocentric tests"""

import unittest

from geographiclib.geocentric import Geocentric


class GeocentricTest(unittest.TestCase):
  """Geocentric test suite"""

  def setUp(self):
      """Instantiate WGS84 Geocentric instance"""
      self.geoc = Geocentric.WGS84()

  def test_forward_origin(self):
    """Test Forward at equator/prime meridian"""
    result = self.geoc.Forward(0, 0, 0)
    self.assertAlmostEqual(result['X'], 6378137.0, places=3)
    self.assertAlmostEqual(result['Y'], 0.0, places=10)
    self.assertAlmostEqual(result['Z'], 0.0, places=10)

  def test_forward_north_pole(self):
    """Test Forward at north pole"""
    result = self.geoc.Forward(90, 0, 0)
    self.assertAlmostEqual(result['X'], 0.0, places=3)
    self.assertAlmostEqual(result['Y'], 0.0, places=3)
    # Polar radius b = a * (1 - f)
    b = self.geoc.a * (1 - self.geoc.f)
    self.assertAlmostEqual(result['Z'], b, places=3)

  def test_forward_south_pole(self):
    """Test Forward at south pole"""
    result = self.geoc.Forward(-90, 0, 0)
    self.assertAlmostEqual(result['X'], 0.0, places=3)
    self.assertAlmostEqual(result['Y'], 0.0, places=3)
    b = self.geoc.a * (1 - self.geoc.f)
    self.assertAlmostEqual(result['Z'], -b, places=3)

  def test_forward_with_height(self):
    """Test Forward with non-zero height"""
    h = 1000.0  # 1 km altitude
    result = self.geoc.Forward(0, 0, h)
    self.assertAlmostEqual(result['X'], 6378137.0 + h, places=3)
    self.assertAlmostEqual(result['Y'], 0.0, places=10)
    self.assertAlmostEqual(result['Z'], 0.0, places=10)

  def test_roundtrip_equator(self):
    """Test Forward->Reverse roundtrip at equator"""
    lat, lon, h = 0.0, 45.0, 100.0
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_roundtrip_midlatitude(self):
    """Test Forward->Reverse roundtrip at mid-latitude"""
    lat, lon, h = 48.8566, 2.3522, 35.0  # Paris
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_roundtrip_high_latitude(self):
    """Test Forward->Reverse roundtrip at high latitude"""
    lat, lon, h = 78.2232, 15.6267, 10.0  # Longyearbyen, Svalbard
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_roundtrip_negative_longitude(self):
    """Test Forward->Reverse roundtrip with negative longitude"""
    lat, lon, h = 40.7128, -74.0060, 10.0  # New York
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_roundtrip_southern_hemisphere(self):
    """Test Forward->Reverse roundtrip in southern hemisphere"""
    lat, lon, h = -33.8688, 151.2093, 58.0  # Sydney
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_roundtrip_high_altitude(self):
    """Test Forward->Reverse roundtrip at high altitude"""
    lat, lon, h = 0.0, 0.0, 35786000.0  # Geostationary orbit
    fwd = self.geoc.Forward(lat, lon, h)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=3)

  def test_roundtrip_poles(self):
    """Test Forward->Reverse roundtrip at poles"""
    for lat in [90.0, -90.0]:
      h = 0.0
      fwd = self.geoc.Forward(lat, 0, h)
      rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
      self.assertAlmostEqual(rev['lat'], lat, places=10)
      # Longitude is undefined at poles, so we don't check it
      self.assertAlmostEqual(rev['h'], h, places=6)

  def test_reverse_origin(self):
    """Test Reverse at earth center returns pole"""
    rev = self.geoc.Reverse(0, 0, 0)
    # At origin, lat should be 90 (north pole), h negative
    self.assertAlmostEqual(rev['lat'], 90.0, places=10)
    self.assertAlmostEqual(rev['lon'], 0.0, places=10)

  def test_spherical_ellipsoid(self):
    """Test with spherical earth (f=0)"""
    R = 6371000.0  # Mean earth radius
    sphere = Geocentric(R, 0)

    # Forward at equator
    fwd = sphere.Forward(0, 0, 0)
    self.assertAlmostEqual(fwd['X'], R, places=3)
    self.assertAlmostEqual(fwd['Y'], 0.0, places=10)
    self.assertAlmostEqual(fwd['Z'], 0.0, places=10)

    # Roundtrip
    lat, lon, h = 45.0, 90.0, 1000.0
    fwd = sphere.Forward(lat, lon, h)
    rev = sphere.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_custom_ellipsoid(self):
    """Test with custom ellipsoid"""
    # GRS80 ellipsoid
    a = 6378137.0
    f = 1/298.257222101
    grs80 = Geocentric(a, f)

    lat, lon, h = 45.0, 45.0, 500.0
    fwd = grs80.Forward(lat, lon, h)
    rev = grs80.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
    self.assertAlmostEqual(rev['lat'], lat, places=10)
    self.assertAlmostEqual(rev['lon'], lon, places=10)
    self.assertAlmostEqual(rev['h'], h, places=6)

  def test_invalid_radius(self):
    """Test that invalid radius raises ValueError"""
    with self.assertRaises(ValueError):
      Geocentric(-1, 0.003)
    with self.assertRaises(ValueError):
      Geocentric(0, 0.003)

  def test_invalid_flattening(self):
    """Test that invalid flattening raises ValueError"""
    # f >= 1 means polar semi-axis is not positive
    with self.assertRaises(ValueError):
      Geocentric(6378137, 1.0)
    with self.assertRaises(ValueError):
      Geocentric(6378137, 2.0)

  def test_dateline_crossing(self):
    """Test coordinates near dateline"""
    for lon in [179.9, 180.0, -179.9, -180.0]:
      lat, h = 0.0, 0.0
      fwd = self.geoc.Forward(lat, lon, h)
      rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'])
      self.assertAlmostEqual(rev['lat'], lat, places=10)
      # Normalize longitude for comparison
      expected_lon = lon if -180 <= lon <= 180 else lon - 360 if lon > 180 else lon + 360
      self.assertAlmostEqual(rev['lon'], expected_lon, places=10)
      self.assertAlmostEqual(rev['h'], h, places=6)

  def test_equatorial_radius(self):
    """Test EquatorialRadius inspector method"""
    self.assertEqual(self.geoc.EquatorialRadius(), self.geoc.a)

  def test_flattening(self):
    """Test Flattening inspector method"""
    self.assertEqual(self.geoc.Flattening(), self.geoc.f)

  def test_forward_with_rotation_matrix(self):
    """Test Forward with M=True returns rotation matrix"""
    result = self.geoc.Forward(45.0, 90.0, 0, M=True)
    self.assertIn('M', result)
    self.assertEqual(len(result['M']), 9)

  def test_reverse_with_rotation_matrix(self):
    """Test Reverse with M=True returns rotation matrix"""
    result = self.geoc.Reverse(6378137.0, 0, 0, M=True)
    self.assertIn('M', result)
    self.assertEqual(len(result['M']), 9)

  def test_forward_without_m_no_matrix(self):
    """Test Forward without M parameter does not include rotation matrix"""
    result = self.geoc.Forward(45.0, 90.0, 0)
    self.assertNotIn('M', result)

  def test_reverse_without_m_no_matrix(self):
    """Test Reverse without M parameter does not include rotation matrix"""
    result = self.geoc.Reverse(6378137.0, 0, 0)
    self.assertNotIn('M', result)

  def test_rotation_matrix_orthonormal(self):
    """Test that rotation matrix is orthonormal (M * M^T = I)"""
    import math
    result = self.geoc.Forward(45.0, 45.0, 0, M=True)
    M = result['M']

    # Check that M * M^T = I (within numerical precision)
    # M is row-major: M[0:3] is row 0, M[3:6] is row 1, M[6:9] is row 2
    for i in range(3):
      for j in range(3):
        # Compute (M * M^T)[i,j] = sum_k M[i,k] * M[j,k]
        dot = sum(M[i*3 + k] * M[j*3 + k] for k in range(3))
        expected = 1.0 if i == j else 0.0
        self.assertAlmostEqual(dot, expected, places=14,
          msg=f"M*M^T[{i},{j}] = {dot}, expected {expected}")

  def test_rotate_unrotate_inverse(self):
    """Test that _Rotate and _Unrotate are inverses"""
    result = self.geoc.Forward(30.0, 60.0, 1000.0, M=True)
    M = result['M']

    # Apply rotation then inverse
    x, y, z = 100.0, 200.0, 300.0
    X, Y, Z = Geocentric._Rotate(M, x, y, z)
    x2, y2, z2 = Geocentric._Unrotate(M, X, Y, Z)

    self.assertAlmostEqual(x2, x, places=10)
    self.assertAlmostEqual(y2, y, places=10)
    self.assertAlmostEqual(z2, z, places=10)

  def test_rotation_equator_prime_meridian(self):
    """Test rotation matrix at equator/prime meridian"""
    result = self.geoc.Forward(0, 0, 0, M=True)
    M = result['M']
    # At lat=0, lon=0:
    # East (local X) -> -Y geocentric (slam=0, clam=1, so M[0]=-0, M[1]=1, M[2]=0)
    # North (local Y) -> -Z geocentric direction for north
    # Up (local Z) -> +X geocentric
    # M[0], M[1], M[2] = -slam, clam, 0 = 0, 1, 0
    self.assertAlmostEqual(M[0], 0.0, places=14)  # -slam
    self.assertAlmostEqual(M[1], 1.0, places=14)  # clam
    self.assertAlmostEqual(M[2], 0.0, places=14)  # 0

  def test_rotation_north_pole(self):
    """Test rotation matrix at north pole"""
    result = self.geoc.Forward(90, 0, 0, M=True)
    M = result['M']
    # At north pole: sphi=1, cphi=0, slam=0, clam=1
    # Up (local Z) should point to +Z geocentric
    self.assertAlmostEqual(M[8], 1.0, places=14)  # sphi

  def test_forward_reverse_rotation_matrices_match(self):
    """Test that Forward and Reverse produce the same rotation matrix"""
    lat, lon, h = 45.0, 90.0, 1000.0
    fwd = self.geoc.Forward(lat, lon, h, M=True)
    rev = self.geoc.Reverse(fwd['X'], fwd['Y'], fwd['Z'], M=True)

    for i in range(9):
      self.assertAlmostEqual(fwd['M'][i], rev['M'][i], places=10,
        msg=f"M[{i}] mismatch: Forward={fwd['M'][i]}, Reverse={rev['M'][i]}")
