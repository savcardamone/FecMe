from FecMe.LDPC import LDPC

class NRLDPC(LDPC):
    """New Radio LDPC Encode/Decode.
    """
    
    def __init__(self, BGN=1):
        """Class constructor.
        """
        self.BGN = BGN
        
    @property
    def BGN(self):
        """Getter for the base graph number.
        """
        return self.__BGN

    @BGN.setter
    def BGN(self, BGN):
        """Setter for the base graph number.
        """
        if BGN == 1 or BGN == 2:
            self.__BGN = BGN
        else:
            raise ValueError("Base graph number {0} is not supported.".format(BGN))
            
    @property
    def Kcb(self):
        """Return the maximum code block size.
        See 38.212 Section 5.2.2
        """
        if self.BGN == 1:
            return 8448
        elif self.BGN == 2:
            return 3840
        else:
            raise ValueError("Can't associate a code block size with base graph number {0}".format(self.BGN))
    
    

    def encode(self, b):
        """Take a bitstring and encode it. We return the fully rate-matched output, g, where
        all segmented code blocks are concatenated.
        See 38.212 Section 5.5 for the end point of this function.
        """

        B = len(b)
