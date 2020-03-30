import os.path
import unittest
import h5py
import numpy as np
from scampy import reads

DATA_DIR = os.path.join(os.path.expanduser("~"), "scampy_test_data")


class TestCITSReads(unittest.TestCase):
    """ Tests for reads. Needs scampy_test_data in HOME folder """

    def assertCloseToBenchmark(self, data, data_type):
        with h5py.File(os.path.join(DATA_DIR, "benchmarks.h5"), "r") as h5file:
            benchmark = h5file[data_type]
            np.testing.assert_allclose(data, benchmark)

    def test_reading_topo(self):
        """ Tests that reading a Topo.txt file works """
        file_path = os.path.join(DATA_DIR, "Topo.txt")
        topo = reads.readTopo(file_path)
        self.assertCloseToBenchmark(topo, "Topo")

    def test_reading_ASCII(self):
        """ Tests that reading an ASCII CITS works. """
        file_path = os.path.join(DATA_DIR, "CITS.asc")
        topo, m_data, channelList, m_params = reads.readCitsAscii(file_path)
        self.assertCloseToBenchmark(m_data, "asc")

    def test_reading_3ds(self):
        """ Tests that reading a 3ds CITS works. """
        file_path = os.path.join(DATA_DIR, "CITS.3ds")
        topo, m_data, channelList, m_params = reads.readCits3dsBin(
            file_path, zSpectro=False
        )
        self.assertCloseToBenchmark(m_data, "3ds")

    def test_reading_sm4(self):
        """ Tests that reading a SM4 CITS works. """
        file_path = os.path.join(DATA_DIR, "CITS.sm4")
        topo, m_data, channelList, m_params, average = reads.readCitsSm4Bin(file_path)
        self.assertCloseToBenchmark(m_data, "sm4")


if __name__ == "__main__":
    unittest.main()
