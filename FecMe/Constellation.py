from math import sqrt, ceil
from numpy import zeros, ones

def ConstellationFactory(constellation_type):
    """Constellation factory object creation. Create a constellation object of the
    user-specified type.
    """
    
    if constellation_type == "QPSK":
        return QPSK(constellation_type)
    else:
        raise ValueError("Unsupported constellation type: {0}".format(constellation_type))

    
class Constellation(object):
    """Base class from which constellations will be derived. Effectively abstract
    with the exception of some functionality for the sake of informative reporting.
    """
    
    def __init__(self, constellation_type):
        self.name = constellation_type
        
    def __str__(self):
        return "Constellation object for {0}".format(self.name)

    def map(self, bitstring):
        raise NotImplementedError("Mapping for {0} constellation is unsupported.".format(self.name))

    def demap(self, constellations):
        raise NotImplementedError("Demapping for {0} constellation is unsupported.".format(self.name))

class QPSK(Constellation):
    """QPSK constellation class.
    """
    
    def __init__(self, constellation_type):
        super(QPSK, self).__init__(constellation_type)
        
    def map(self, bitstring):
        """Map a bitstring to a series of points in a QPSK constellation.
        """

        num_constellations = int(ceil(len(bitstring) / 2))        
        constellations = zeros((num_constellations,1), dtype=complex)
        
        for i in range(num_constellations):
            # If we have an odd number of bits to map then we can just handle the index error we
            # get from accessing past the end of the array by setting the imaginary part of the
            # constellation to zero
            try:
                constellations[i] = ((2*bitstring[2*i] - 1) + 1j*(2*bitstring[2*i+1] - 1)) / sqrt(2)
            except IndexError:
                constellations[i] = (2*bitstring[2*i] - 1) / sqrt(2)
                
        return constellations

    def demap(self, constellations):
        """Demap constellation points to LLRs for underlying bits.
        """
    

if __name__ == "__main__":

    my_qpsk = ConstellationFactory("QPSK")
    print(my_qpsk)

    my_bitstring = ones((19,1), dtype=int)
    my_constellations = my_qpsk.map(my_bitstring)
    print(my_constellations)
    my_recovered = my_qpsk.demap(my_constellations)
