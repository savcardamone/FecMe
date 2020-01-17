from os.path import abspath, dirname
from math import inf
from numpy import array, zeros, ceil, nonzero, min, where, argwhere, concatenate, load, identity, roll

from FecMe.CRC import polynomials, checksum, check

class NRLDPC():
    """New Radio LDPC Encode/Decode.
    """

    lifting_sizes = array(
        [[2,4,8,16,32,64,128,256],
         [3,6,12,24,48,96,192,384],
         [5,10,20,40,80,160,320,0],
         [7,14,28,56,112,224,0,0],
         [9,18,36,72,144,288,0,0],
         [11,22,44,88,176,352,0,0],
         [13,26,52,104,208,0,0,0],
         [15,30,60,120,240,0,0,0]], dtype=int)

    def __init__(self, A, BGN=1):
        """Class constructor.
        """
        self.A = A
        self.BGN = BGN
        self.BG_cache = None
        self.PCM_cache = None
        
    def __str__(self):
        """String representation of NRLDPC object.
        """
        return "NRLDPC: {0} input bits, Base-Graph Number {1}".format(self.A, self.BGN)
        
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
        if A > 0:
            self.__A = A
        else:
            raise ValueError("Number of bits for coding must be greater than zero.")
        
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
        else:
            return 3840

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
            return 1
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
        return self.Bprime // self.C

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

    @property
    def Zc(self):
        """Return the lifting size for the base graph.
        """
        local_lifting_sizes = NRLDPC.lifting_sizes
        # Pull out the indices whose entries don't satisfy the inequality Kb.Zc >= K'
        sub_threshold_indices = local_lifting_sizes < self.Kprime / self.Kb
        # Set subthreshold elements to zero
        local_lifting_sizes[sub_threshold_indices] = 0
        # Now pull out the smallest non-zero element from the lifting sizes and that's our
        # desired lifting size
        lifting_set, z = where(local_lifting_sizes == min(local_lifting_sizes[nonzero(local_lifting_sizes)]))

        return int(local_lifting_sizes[lifting_set,z])
        
    @property
    def LiftingSet(self):
        """Return the lifting set for the base graph.
        """
        local_lifting_sizes = NRLDPC.lifting_sizes
        lifting_set = argwhere(local_lifting_sizes == self.Zc)[0][0]

        return int(lifting_set)
    
    @property
    def K(self):
        """Return the number of bits the encoded codeblock will comprise.
        """
        if self.BGN == 1:
            return 22 * self.Zc
        else:
            return 10 * self.Zc

    @property
    def BG(self):
        """Return the base graph from which the LDPC PCM will be constructed. If the base graph
        hasn't yet been formed, then explicitly do so, otherwise return the cached base graph.
        """
        if self.BG_cache is None:
            # Get the path to FecMe/FecMe so we can find the npz file
            path = abspath(dirname(__file__))
            # Select the appropriate base graph
            BGs = load('{0}/NRLDPC_Base_Graphs.npz'.format(path))
            self.BG_cache = BGs['BaseGraph{0}_LiftingSet{1}'.format(self.BGN, self.LiftingSet)]

        return self.BG_cache
        
    @property
    def PCM(self):
        """Return the parity check matrix (PCM) for the LDPC code. If it hasn't yet been created,
        explicitly make it, otherwise return the cached PCM.
        See 38.212 Section 5.3.2, Item (3) for details about construction
        """
        # If we don't have a PCM in cache, then explicitly construct it
        if self.PCM_cache is None:
            
            bg_nrows = self.BG.shape[0]
            bg_ncols = self.BG.shape[1]

            # Now our PCM cache contains something, so we won't enter this branch again
            self.PCM_cache = zeros((bg_nrows*self.Zc,bg_ncols*self.Zc), dtype=int)

            # Loop over each entry in base graph and 
            for irow in range(bg_nrows):
                for icol in range(bg_ncols):

                    # If the base graph entry is zero, the PCM receives a null matrix
                    if self.BG[irow,icol] == 0:
                        self.PCM_cache[irow*self.Zc:(irow+1)*self.Zc,icol*self.Zc:(icol+1)*self.Zc] = zeros((self.Zc,self.Zc), dtype=int)
                    # If the base graph entry is nonzero, the PCM receives a cyclically shifted
                    # identity matrix
                    else:
                        shift = self.BG[irow,icol] % self.Zc
                        self.PCM_cache[irow*self.Zc:(irow+1)*self.Zc,icol*self.Zc:(icol+1)*self.Zc] = roll(identity(self.Zc, dtype=int), shift, axis=1)

        return self.PCM_cache
            
        
    ### ====================================================================================
    ###                                     Methods
    ### ====================================================================================

    def segmentation(self, b):
        """Segment the transport block (with attached CRC bits) into codeblocks which
        are to be encoded separately.
        See 38.212 Section 5.2.2. for details.
        """
        # Initialise our segmented codeblocks with zeros
        c = zeros((self.C,self.K), dtype=int)

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
                p = checksum(c[r,0:(self.Kprime - self.L)], polynomial='CRC24B', checksum_fill=0)
                # Append the checksum bits to the message bits in the codeblock
                for k in range((self.Kprime-self.L), self.Kprime):
                    c[r,k] = p[k + self.L - self.Kprime]

            # Pad the remaining bits in the codeblock with filler bits
            for k in range(self.Kprime, self.K):
                c[r,k] = -1

        return c
    
    def parity(self, c, verify=False):
        """Encode each codeblock with the lifted base graph and return the systematic 
        codewords.
        See 38.212 Section 5.3.2. for details, although the means by which the parity bits are
        derived is not discussed.
        """

        # Verify that the produced codeword spans the nullspace of the PCM if we're explicitly
        # verifying the encoding
        if verify == True:
            pass
        
        pass


    def encode(self, a):
        """Take a bitstring and encode it. We return the fully rate-matched output, g, where
        all segmented code blocks are concatenated.
        See 38.212 Section 5.5 for the end point of this function.
        """


        def rate(self, d):
            def select(self):
                pass
            def interleave(self):
                pass
            pass

        def concatenation(self, f):
            pass
        
        # Make sure we're not trying to encode something that doesn't belong
        if len(a) != self.A:
            raise ValueError("Encoder has been parameterised for {0} bits, but {1} have been passed".format(self.A, len(a)))

        # QQ: Not sure if this is the right polynomial. Can't find it in the spec...
        b = concatenate((a, checksum(a, polynomial='CRC24A', checksum_fill=0)))
        # Transport block segmentation
        c = self.segmentation(b)
        # Generate parity bits for each codeblock
        d = parity(self, c)
        # Rate matching
        f = rate(self, d)
        # Codeblock concatenation
        g = concatenation(self, f)

        return g
