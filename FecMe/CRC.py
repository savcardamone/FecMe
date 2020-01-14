### FILE: CRC.py
### AUTHOR: Salvatore Cardamone
### DESCRIPTION: Various cyclic-redundancy check (CRC) functionalities and polynomials

from numpy import array, zeros, ones, trim_zeros, concatenate, argmax

# Generator polynomials taken from Section 5.1 of 38.212
# Leading element in array corresponds to highest power term, final element is the zeroth power
polynomials = {
    'CRC24A' : array([1,1,0,0,0,0,1,1,0,0,1,0,0,1,1,0,0,1,1,1,1,1,0,1,1], dtype=int),
    'CRC24B' : array([1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,1,1], dtype=int),
    'CRC24C' : array([1,1,0,1,1,0,0,1,0,1,0,1,1,0,0,0,1,0,0,0,1,0,1,1,1], dtype=int),
    'CRC16'  : array([1,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1], dtype=int),
    'CRC11'  : array([1,1,1,0,0,0,1,0,0,0,0,1], dtype=int),
    'CRC6'   : array([1,1,0,0,0,0,1], dtype=int)
}

def checksum(a, polynomial='CRC24A', checksum_fill=0):
    """Take a bitstring and compute the CRC checksum to be appended to the bitstring.
    Return CRC checksum -- the calling function can append them if needs be.
    """

    # We have no need to consider the leading zeros in the input bitstring so
    # just hack them off right at the start
    #a = trim_zeros(a, trim='f')

    A = len(a)
    g = polynomials.get(polynomial)
    L = len(g)

    # Add overflow portion to work array where the checksum will end up
    if checksum_fill == 0:
        work = concatenate((a, zeros((L-1,), dtype=int)))
    elif checksum_fill == 1:
        work = concatenate((a, ones((L-1,),  dtype=int)))
    else:
        raise ValueError("Can't initialise the checksum with {0}s".format(checksum_fill))

    # Continue until all of the message part of the work array is zero
    while any(work[0:A]) == 1:

        # Find the leading non-zero element in the work array -- that's where we apply the
        # generator polynomial
        current_shift = argmax(work > 0)

        # Loop over bits in the generator polynomial and xor with bit directly above
        # in the work array
        work[range(current_shift,current_shift+L)] ^= g[:]

    # Return the CRC checksum
    return work[A:]

def check(b, polynomial='CRC24A', checksum_fill=0):
    """Verify whether the received CRC checksum is valid given the generating polynomial.
    """
    
    L = len(polynomials.get(polynomial))
    A = len(b) - L

    # Compute the checksum from the received message bits
    crc_checksum = checksum(b[0:A], polynomial=polynomial, checksum_fill=checksum_fill)
    # Pull out what we received as the checksum
    rx_crc_checksum = b[A:]

    # Verify whether or not the checksums match and flag to the caller
    return False if any(crc_checksum ^ rx_crc_checksum) == ~checksum_fill else True
    
if __name__ == "__main__":

    test_bitstring = array([0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,1,1,1,1,1,1], dtype=int)
    crc_checksum = checksum(test_bitstring, polynomial="CRC16")
    print("CRC Checksum: {0}".format(crc_checksum))
    crc_pass = check(concatenate((test_bitstring, crc_checksum)), polynomial="CRC16")
    print("CRC Pass: {0}".format(crc_pass))
