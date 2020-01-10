from numpy import array
from scipy.linalg import lu

class LDPC():

    def __init__(self, pcm):
        self.pcm = pcm

    def __str__(self):
        return "({0} x {1}) PCM Matrix\n{2}".format(self.n - self.k, self.n, self.pcm)
        
    ### ====================================================================================
    ###                                     Parameters
    ### ====================================================================================

    @property
    def pcm(self):
        """Getter for the parity check matrix defining the LDPC code.
        """
        return self.__pcm

    @pcm.setter
    def pcm(self, pcm):
        """Setter for the parity check matrix defining the LDPC code. Also preprocess the
        PCM to get it into row-reduced form.
        """
        pl, u = lu(pcm, permute_l=True)
        self.__pcm = (u % 2).astype(int)
    
    ### ====================================================================================
    ###                                 Dependent Variables
    ### ====================================================================================

    @property
    def n(self):
        """Return the length of the LDPC code.
        """
        return self.pcm.shape[1]

    @property
    def k(self):
        """Return the dimension of the LDPC code.
        """
        return self.n - self.pcm.shape[0]
    
    ### ====================================================================================
    ###                                     Methods
    ### ====================================================================================

    def encode():
        pass
    

if __name__ == "__main__":

    pcm = array([[1,1,1,1,0,0],[0,0,1,1,0,1],[1,0,0,1,1,0]], dtype=bool)
    ldpc = LDPC(pcm)
    print(ldpc)
