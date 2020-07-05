"""
Diffuse information across networks.
"""

import numpy as np
import pandas as pd

from scipy.sparse.linalg import lgmres
from scipy import sparse

from kinapp_helper import InputValidator

class Diffusion:
    """
    Propagate label across a graph.

    Attributes
    ----------
    attr1 : type
        description
    """

    def __init__(self, network, labels):
        """
        Inits diffusion experiment with network and starting labels.

        Parameters
        ----------
        network : sparse coo matrix
            graph represented as an adjacency matrix
        labels : list
            indeces of nodes that carry the initial label
        """

        self.network = network
        self.labels = labels

    def get_label_indices(self, labels):
        """
        Converts protein ids to their network indices.

        Returns
        -------
        labels_indexed : list[int]
            position of proteins in the network matrix
        """
        hugo_ids = {k : v for v, k in enumerate(self.network.protein_list)}

        labels_indexed = []
        for label in labels:
            labels_indexed.append(hugo_ids[label])
        return labels_indexed

    def diffuse(self):
        """
        Diffuse labels across network.
        """

        lpp = sparse.csgraph.laplacian(self.network.adj_matrix)
        alpha = 1 / float(np.max(np.sum(np.abs(lpp), axis=0)))
        ident = sparse.csc_matrix(np.eye(self.network.adj_matrix.shape[0]))
        ps = ident + alpha*lpp

        initial_label_state = np.zeros(self.network.adj_matrix.shape[0])
        label_indices = self.get_label_indices(self.labels)
        initial_label_state[label_indices] = 1

        #diffuse
        diff_out = lgmres(ps, initial_label_state, maxiter=1000, atol=0.001)#, atol=0, maxiter=1e3)
        final_label_state = diff_out[0]
        result = DiffusionResult(final_label_state,
                                 initial_label_state,
                                 self.network.protein_list)
        return result

class DiffusionResult:
    """
    Class to namespace functions for dealing with diffusion output.
    """

    def __init__(self, result, labels, proteins):
        """
        Inits instance with diffusion output and corresponding protein names.

        Parameters
        ----------
        result : numpy array
            Output of scipy.sparse.linalg.lgmres call; diffusion result
        labels : list[str]
            List of input labels in the diffusion
        """
        self.initial_labels = labels
        self.final_labels = result
        self.proteins = proteins
        self.protein_to_final_label = dict(zip(proteins, result))


        #print(len(result), len(labels), len(self.proteins))
        #print(result, labels, self.proteins)

    def get_as_pandas_df(self, include_z_score = False):
        """
        Formats diffusion output as pandas df
        """
        self.result = pd.DataFrame({'protein':self.proteins,
                                    'initial_label':self.initial_labels,
                                    'final_label':self.final_labels})
        self.result.sort_values(by='final_label', ascending=False, inplace=True)

        if include_z_score:
            z_score_result = self.result[self.result.initial_label==0].copy()
            z_score_result['z_score'] = (z_score_result.final_label
                                         - z_score_result.final_label.mean())/z_score_result.final_label.std()
            return z_score_result
        return self.result

    def get_result_for_protein(self, protein):
        """
        Get post-diffusion score for specific protein.
        """
        return self.protein_to_final_label[protein]
