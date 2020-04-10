import unittest
import numpy as np
from scampy import processing

YPX = 10
XPX = 10


class TestProcessing(unittest.TestCase):
    """ Tests for processing. """

    def gen_dummy_cits_data(self):
        return np.random.random((2, YPX, XPX, 100))

    def gen_dummy_topo(self):
        return np.random.random((YPX, XPX))

    def test_directionAverageCITS(self):
        """ Tests that the direction average returns the approriate shape. """
        cits_data = self.gen_dummy_cits_data()
        binning = 2
        initial_shape = cits_data.shape
        avg_y = processing.directionAverageCITS(cits_data, binning, "y")
        self.assertEqual(
            avg_y.shape,
            (
                initial_shape[0],
                initial_shape[1] // binning,
                initial_shape[2],
                initial_shape[3],
            ),
        )
        avg_x = processing.directionAverageCITS(cits_data, binning, "x")
        self.assertEqual(
            avg_x.shape,
            (
                initial_shape[0],
                initial_shape[1],
                initial_shape[2] // binning,
                initial_shape[3],
            ),
        )

    def test_normalizeDOS(self):
        """ Tests that normalizeDOS returns the approriate shape. """
        cits_data = self.gen_dummy_cits_data()
        norm_cits = processing.normalizeDOS(cits_data)
        self.assertEqual(norm_cits.shape, cits_data.shape)

    def test_levelTopo(self):
        topo = self.gen_dummy_topo()
        leveled_topo = processing.levelTopo(topo)
        self.assertEqual(leveled_topo.shape, topo.shape)

    def test_extractSlope(self):
        dummy_channel = self.gen_dummy_cits_data()[0]
        dummy_topo = self.gen_dummy_topo()
        slope_data, coef_data, zg = processing.extractSlope(
            dummy_topo, dummy_channel, 0.1, 0.001
        )
        self.assertEqual(slope_data.shape, dummy_channel.shape)
        self.assertEqual(coef_data.shape, dummy_channel.shape)
        self.assertEqual(zg.shape, dummy_channel.shape)

    def test_findPixelsOnLine_45deg(self):
        """ The algo should find pixels of coordinates (i,i) for a line at 45 degrees."""
        n_pixels = 5
        xi, xf = 0, n_pixels
        yi, yf = 0, n_pixels
        pixels_x, pixels_y = processing.findPixelsOnLine(xi, xf, yi, yf)
        for i in range(n_pixels):
            self.assertEqual(i, pixels_x[i])
            self.assertEqual(i, pixels_y[i])

    def test_findPixelsOnLine_vertical(self):
        """ The algo should find pixels of coordinates (0,i) for a vertical line."""
        n_pixels = 5
        xi, xf = 0, 0
        yi, yf = 0, n_pixels
        pixels_x, pixels_y = processing.findPixelsOnLine(xi, xf, yi, yf)
        for i in range(n_pixels):
            self.assertEqual(0, pixels_x[i])
            self.assertEqual(i, pixels_y[i])

    def test_findPixelsOnLine_horizontal(self):
        """ The algo should find pixels of coordinates (i,0) for an horizontal line."""
        n_pixels = 5
        xi, xf = 0, n_pixels
        yi, yf = 0, 0
        pixels_x, pixels_y = processing.findPixelsOnLine(xi, xf, yi, yf)
        for i in range(n_pixels):
            self.assertEqual(i, pixels_x[i])
            self.assertEqual(0, pixels_y[i])

    def test_findPixelsOnLine_any(self):
        """ The algo should a least include the first and the last pixel for any line. """
        n_pixels = 5
        xi, xf = np.random.randint(0, n_pixels, size=2)
        yi, yf = np.random.randint(0, n_pixels, size=2)
        pixels_x, pixels_y = processing.findPixelsOnLine(xi, xf, yi, yf)
        self.assertEqual(xi, pixels_x[0])
        self.assertEqual(xf, pixels_x[-1])
        self.assertEqual(yi, pixels_y[0])
        self.assertEqual(yf, pixels_y[-1])

    def test_stringify(self):
        """ Tests that an array is correctly converted in string. """
        test_string = "Array"
        array_to_str = np.array([ord(x) for x in test_string], dtype=int)
        self.assertEqual(processing.stringify(array_to_str), test_string)


if __name__ == "__main__":
    unittest.main()
