### FILE: test_NRLDPC.py
### AUTHOR: Salvatore Cardamone
### DESCRIPTION: Verify our NR LDPC implementation against golden data

import unittest
from numpy import load

from FecMe.CRC import checksum
from FecMe.NRLDPC import NRLDPC

class TestNRLDPC(unittest.TestCase):
    """NR LDPC Unit testing.
    """
    

    def __init__(self, *args, **kwargs):
        """Class constructor. Load all the test vectors once.
        """
        super(TestNRLDPC, self).__init__(*args, **kwargs)
        self.test_vectors = load('test/test_NRLDPC.npz')
    
    def test_1(self):
        """Test the transport block CRC against golden data.
        """
        
        crc_bits = checksum(self.test_vectors['data_in'], polynomial='CRC24A', checksum_fill=0)
        golden_crc_bits = self.test_vectors['data_in_crc'][-24:]

    def test_2(self):
        """Test the NR LDPC encoder.
        """
        pass
