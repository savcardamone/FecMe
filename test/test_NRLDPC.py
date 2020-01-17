### FILE: test_NRLDPC.py
### AUTHOR: Salvatore Cardamone
### DESCRIPTION: Verify our NR LDPC implementation against golden data

import unittest
from numpy import load, array

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
        golden_crc_out = self.test_vectors['crc_out'][-24:]
        crc_out = checksum(self.test_vectors['data_in'], polynomial='CRC24A', checksum_fill=0)
        self.assertEqual(crc_out.tolist(), golden_crc_out.tolist())
        
    def test_2(self):
        """Test the NR LDPC segmentation against golden data.
        """
        A = len(self.test_vectors['data_in'])
        ldpc = NRLDPC(A)

        # Had to dump each codeblock into a column in MATLAB, so the codeblocks are interleaved in
        # the numpy array. Just deinterleave them and get them into the correct format
        golden_seg_out = array([self.test_vectors['seg_out'][::2], self.test_vectors['seg_out'][1::2]])
        seg_out = ldpc.segmentation(self.test_vectors['crc_out'])
        test_var = self.assertEqual(seg_out.tolist(), golden_seg_out.tolist())
        
    def test_3(self):
        """Test the NR LDPC codeword creation against golden data.
        """
        pass

    def test_4(self):
        """Test the full NR LDPC encode chain against golden data.
        """
        pass
