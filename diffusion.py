from scipy import sparse
import numpy as np
from scipy.sparse.linalg import lgmres


class Diffusion:
    """
    Propagate label across a graph
    
    Attributes
    ----------
    attr1 : type
        description

    Methods
    ------

    """

    def __init__(self, network, labels):
        """
        Parameters
        ---------
        network : sparse coo matrix
            graph represented as an adjacency matrix
        labels : list
            indeces of nodes that carry the initial label
        """
        
        self.network = network
        self.labels = labels

    def diffuse(self):
        """
        Diffuse labels across network
        """
        
        lpp = self.get_lpp()
        alpha = self.get_diffusion_parameter(lpp)
        ps = self.get_ps(alpha, lpp) 
        
        y = np.zeros([self.network.shape[0],1])
        y[0] = 1
        #diffuse
        result = lgmres(ps, y, maxiter=1000, atol=0.001)#, atol=0, maxiter=1e3)
        return result
    
    def get_ps(self, lpp, alpha):

        ident = sparse.csc_matrix(np.eye(self.network.shape[0]))
        ps = ident + alpha*lpp
        return ps

    def get_lpp(self):
        """
        Calculates the network Laplacian
        """
        return sparse.csgraph.laplacian(self.network)

    def get_diffusion_parameter(self, lpp):
        diffusionParameter = 1 / float( np.max( np.sum(np.abs(lpp), axis=0) ) )
        return diffusionParameter

        

