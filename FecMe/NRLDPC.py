from numpy import ceil

from FecMe.CRC import polynomials, checksum
from FecMe.LDPC import LDPC

class NRLDPC(LDPC):
    """New Radio LDPC Encode/Decode.
    """
    
    def __init__(self, A, BGN=1):
        """Class constructor.
        """
        self.A = A
        self.BGN = BGN

    ### ====================================================================================
    ###                                     Parameters
    ### ====================================================================================

    @property
    def A(self):
        """Getter for the number of bits in the transport block.
        """
        return self.__A

    @A.setter
    def A(self, A):
        """Setter for the number of bits in the transport block.
        """
        self.__A = A
        
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

    ### ====================================================================================
    ###                                 Dependent Variables
    ### ====================================================================================
        
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

    @property
    def B(self):
        """Return the number of bits in the transport block plus the appended CRC.
        """
        return self.A + 24

    @property
    def L(self):
        """Return the number of bits that get appended to each segmented codeblock.
        """
        if self.B <= self.Kcb:
            return 0
        else:
            return 24
        
    @property
    def C(self):
        """Return the number of codeblocks the transport block is to be segmented into.
        """
        if self.B <= self.Kcb:
            return 0
        else:
            return int(ceil(self.B / (self.Kcb - self.L)))

    @property
    def Bprime(self):
        """Return the total number of bits of all codeblocks combined.
        """
        return self.B + self.C * self.L

    @property
    def Kprime(self):
        """Return the number of bits in each codeblock.
        """
        return self.Bprime / self.C

    @property
    def Kb(self):
        """Return Kb, a parameter which helps us select the lifting size.
        """
        if self.BGN == 1:
            return 22
        elif self.BGN == 2:
            if self.B > 640:
                return 10
            elif self.B > 560:
                return 9
            elif self.B > 192:
                return 8
            else:
                return 6
    

    ### ====================================================================================
    ###                                     Methods
    ### ====================================================================================

    def encode(self, a):
        """Take a bitstring and encode it. We return the fully rate-matched output, g, where
        all segmented code blocks are concatenated.
        See 38.212 Section 5.5 for the end point of this function.
        """

        def segmentation(self, b):
            """Segment the transport block (with attached CRC bits) into codeblocks which
            are to be encoded separately.
            See 38.212 Section 5.2.2. for details.
            """
    
            # Initialise our segmented codeblocks with zeros
            c = zeros((self.C,self.Kprime), dtype=int)

            # Counter we manually increment when filling message bits
            s = 0
            # Loop over codeblocks we're segmenting the transport block into
            for r in range(self.C):
                
                # Fill the message bits in the codeblock
                for k in range(self.Kprime - self.L):
                    c[r,k] = b[s]
                    s += 1

                # If we're segmenting the transport block, each codeblock gets its own CRC checksum 
                if self.C > 1:
                    p = CRC.checksum(c[r,self.Kprime - self.L], polynomial='CRC24B', checksum_fill=0)
                    # Append the checksum bits to the message bits in the codeblock
                    for k in range((self.Kprime - self.L):self.Kprime):
                        c[r,k] = p[k + self.L - self.Kprime]

                # Pad the remaining bits in the codeblock with filler bits
                for k in range(self.Kprime:self.K):
                    c[r,k] = -1

            return c

        def parity(self, c):
            pass
        
        # Make sure we're not trying to encode something that doesn't belong
        if len(a) != self.A:
            raise ValueError("Encoder has been parameterised for {0} bits, but {1} have been passed".format(self.A, len(a)))

        # QQ: Not sure if this is the right polynomial. Can't find it in the spec...
        b = concatenate(a, CRC.checksum(a, polynomial='CRC24A', checksum_fill=0))
        # Transport block segmentation
        c = segmentation(self, b)
        # Generate parity bits for each codeblock
        d = parity(self, c)
