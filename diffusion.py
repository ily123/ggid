"""Diffuse information across network.

Usage:

    diffusion_experiment = Diffusion(network, input_nodes)
    result = diffusion.diffuse()
    print(result.head())
"""


import numpy as np
import pandas as pd
from scipy import sparse, stats
from scipy.sparse.linalg import lgmres

pd.options.mode.chained_assignment = None


class Diffusion:
    """Propagate information across a graph."""

    def __init__(self, network, input_nodes):
        """Inits diffusion experiment with network and starting input nodes.

        Parameters
        ----------
        network : sparse coo matrix
            graph represented as an adjacency matrix
        input_nodes: list
            string ids of nodes that carry the initial label
        """
        self.network = network
        self.input_nodes = input_nodes  # this is where information is diffused from

    def get_node_indices(self, proteins):
        """Gets network position (matrix indices) for set of protein ids."""
        hugo_ids = {k: v for v, k in enumerate(self.network.proteins)}
        indices = [hugo_ids[protein] for protein in proteins]
        return indices

    def diffuse(self):
        """Diffuses information from input nodes across the graph."""
        lpp = sparse.csgraph.laplacian(self.network.network)
        alpha = 1 / float(np.max(np.sum(np.abs(lpp), axis=0)))
        ident = sparse.csc_matrix(np.eye(self.network.network.shape[0]))
        ps = ident + alpha * lpp

        initial_state = np.zeros(self.network.network.shape[0])
        input_indices = self.get_node_indices(self.input_nodes)
        initial_state[input_indices] = 1

        diff_out = lgmres(ps, initial_state, maxiter=1000, atol=1e-12)
        final_state = diff_out[0]
        result = DiffusionResult(final_state, initial_state, self.network.proteins)
        return result


class DiffusionResult:
    """Class to namespace functions for dealing with diffusion output."""

    def __init__(self, result, input_nodes, proteins):
        """Inits instance with diffusion output and corresponding protein names.

        Parameters
        ----------
        result : numpy array
            Output of scipy.sparse.linalg.lgmres call; diffusion result
        labels : list[str]
            List of input labels in the diffusion
        """
        self.initial_state = input_nodes
        self.final_state = result
        self.proteins = proteins
        self.protein_final_state = dict(zip(proteins, result))

    def get_result_df(self):
        """Formats diffusion result as pandas df."""
        result = pd.DataFrame(
            {
                "protein": self.proteins,
                "initial_state": self.initial_state,
                "final_state": self.final_state,
            }
        )
        return result

    def get_result_df_with_zscore(self):
        """Formats diffusion result as pandas df and adds zscore & rank column."""
        result = self.get_result_df()
        result = result[result.initial_state == 0]
        result["zscore"] = stats.zscore(result.final_state)
        result["rank"] = result.zscore.rank(ascending=False)
        result.sort_values(by="final_state", ascending=False, inplace=True)
        return result

    def get_result_for_protein(self, protein):
        """Gets post-diffusion score for specific protein."""
        return self.protein_final_state[protein]
