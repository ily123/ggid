print("reloaded")
"""
Calculates GO term similarity between a set of proteins to
produce a similarity matrix, which can then be converted into a
protein-protein network for diffusion.
"""

import multiprocessing
import pickle
import time
from typing import Dict, List, Type, Union

import numpy as np
import pandas as pd
from scipy import sparse

import onto


class SimilarityCalculator2:
    """Calculate Resnik similarity between proteins/entities in a set."""

    def __init__(
        self,
        annotations: onto.Annotations,
        ontology: onto.GoGraph,
        proteins: List[str],
    ) -> None:
        """Inits with annotations, GO graph, and protein list.

        Parameters
        ----------
        annotations : onto.Annotations
            protein annotation corpus
        ontology : onto.GoGraph
            Gene Ontology term graph
        proteins : List[str]
            list of proteins (HUGO ids) to build network for
        """
        self.annotations = annotations.get_as_dict()
        self.ontology = ontology
        self.proteins = proteins
        self.protein_similarity = []
        self._term_similarity = {}

    def filter_proteins(self, min_annotations: int) -> Dict[str, List[str]]:
        """Fileters out under-annotated or missing proteins from the protein list.

        Parameters
        ----------
        min_annotations : int
            min number of annotations a protein must have to be retained in the list

        Returns
        -------
        filtered_out : Dict[str, List]
            dictionary of proteins that got filtered out
        """
        # proteins not in the annotations corpus
        present = self.is_present()
        not_present = set(self.proteins) - set(present)
        # proteins with insufficient annotation count
        underannotated = [
            p for p in self.proteins if len(self.annotations[p]) < min_annotations
        ]
        # proteins to keep
        self.proteins = list(set(present) - set(underannotated))
        filtered_out = dict(not_present=not_present, underannotated=underannotated)
        return filtered_out

    def is_present(self) -> List[str]:
        """Returns list of proteins present in the annotation corpus.

        Returns
        -------
        present : List[str]
            subset of proteins in self.protein that is also found the annotation corpus
        """
        present = list(set(self.annotations) & set(self.proteins))
        return present

    def calculate_similarity(self) -> None:
        """Calculate similarity of all protein pairs in the set."""
        self.proteins = sorted(self.proteins)
        protein_a_index = []
        protein_b_index = []
        similarity = []
        for index_a, protein_a in enumerate(self.proteins):
            for index_b, protein_b in enumerate(self.proteins):
                if index_b >= index_a:
                    sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
                    protein_a_index.append(index_a)
                    protein_b_index.append(index_b)
                    similarity.append(sim)
        similarity = sparse.coo_matrix((similarity, (protein_a_index, protein_b_index)))
        self.protein_similarity = similarity

    def calculate_similarity_two_proteins(
        self, protein_a: str, protein_b: str
    ) -> float:
        """Calculate Resnik (BMA) similarity between two proteins.

        Parameters
        ----------
        protein_a : str
            protein id (HUGO) for first protein
        protein_b : str
            protein id (HUGO) for second protein

        Returns
        -------
        protein_similarity : float
            similarity score, higher is better
        """
        terms_a = self.annotations[protein_a]
        terms_b = self.annotations[protein_b]
        # run go terms through all-vs-all similarity test
        term_similarity = None
        term_sim_vec = []
        for a in terms_a:
            for b in terms_b:
                if (a, b) in self._term_similarity:
                    term_similarity = self._term_similarity[(a, b)]
                else:
                    term_similarity = self.get_term_similarity(a, b)
                    self._term_similarity[(a, b)] = term_similarity
                    self._term_similarity[(b, a)] = term_similarity
                term_sim_vec.append(term_similarity)
        term_sim_matrix = np.array(term_sim_vec).reshape(len(terms_b), len(terms_a))
        # do best match averaging of term similarity scores
        best_match_a = term_sim_matrix.max(axis=0)
        best_match_b = term_sim_matrix.max(axis=1)
        # the BMA score is proxy for protein similarity
        protein_similarity = np.concatenate((best_match_a, best_match_b)).mean()
        return protein_similarity

    def get_term_similarity(self, term1: str, term2: str) -> float:
        """Given two terms, find their similarity.

        Similarity is defined as the specificity of the two terms'
        most specific common ancestor.

        Parameters
        ----------
        term1 : str
            id of the first GO term of interest
        term2 : str
            id of the second GO term of interest

        Returns
        -------
        specificity_mica : float
            most informative common ancestor specificity
        """
        if term1 == term2:
            return self.ontology.nodes[term1].specificity

        ancestors1 = self.ontology.get_full_ancestry(term1)
        ancestors2 = self.ontology.get_full_ancestry(term2)
        shared_ancestors = list(set(ancestors1) & set(ancestors2))
        if len(shared_ancestors) == 0:
            return 0
        if len(shared_ancestors) == 1:
            return self.ontology.nodes[shared_ancestors[0]].specificity
        specificity_mica = 0
        for ancestor in shared_ancestors:
            specificity_mica = max(
                specificity_mica, self.ontology.nodes[ancestor].specificity
            )
        return specificity_mica


class SimilarityMatrix:
    """Similairity matrix and its protein ids."""

    def __init__(self, raw_similarity, proteins):
        """Inits the similarity matrix container.

        Parameters
        ----------
        raw_similarity : numpy 2d array
            Matrix containing protein-protein GO term similarity scores, so that
            raw_similarity[i, j] contains similarity score for proteins i and j.

        proteins : list[str]
            List of proteins in the matrix, sorted so that
            row (column) i in raw_similarity corresponds to protein[i].
        """

        self.raw_similarity = raw_similarity + raw_similarity.T
        self.raw_similarity = self.raw_similarity.todense()
        np.fill_diagonal(self.raw_similarity, 0)
        self.proteins = proteins

    def get_protein_name_by_index(self, protein_index):
        """Get protein name given it's position (index) in network matrix."""

        return self.proteins[protein_index]

    def get_protein_index(self, protein_name):
        """Get protein index (position) in network matrix given protein's name."""

        if protein_name not in self.proteins:
            return None
        else:
            return self.proteins.index(protein_name.upper())

    def get_similarity_vector_for_protein(self, protein_name):
        """Given a protein return it's similarity to other proteins in network."""

        protein_index = self.get_protein_index(protein_name)
        if protein_index:
            similarity_vector = self.raw_similarity[protein_index, :].tolist()[0]
            return pd.DataFrame(
                {"protein": self.proteins, "similarity_score": similarity_vector}
            )
        else:
            return None

    def threshold_matrix(self, n=None):
        """Converts similarity matrix into a network.

        For each protein, sets its N most similar connections to 1,
        and zero-outs the rest.

        Parameters
        ----------
        n : int
            Number of edges to keep for each protein,
            sqrt(network_size) by default.

        Returns
        -------
        network : scipy.sparse.coo.coo_matrix
            An adjacency matrix of 1s and 0s, where 1
            denotes a connection between two protein.
        """

        if not n:
            n = np.ceil(np.sqrt(len(self.proteins)))
        n = int(n)
        adj_matrix = self.raw_similarity.copy()
        net_size = len(self.proteins)
        top_n_edge_cutoff = np.partition(adj_matrix, net_size - n, axis=1)[
            :, net_size - n
        ]
        mask_top_n_edge = adj_matrix >= top_n_edge_cutoff
        adj_matrix[mask_top_n_edge] = 1
        adj_matrix[~mask_top_n_edge] = 0
        self.adj_matrix = adj_matrix

    def enforce_adj_matrix_symmetry(self):
        """Updates the adjacency matrix to be symmetric about the diagonal."""

        symmetric_adj_matrix = self.adj_matrix + self.adj_matrix.T
        symmetric_adj_matrix[symmetric_adj_matrix > 0] = 1
        self.adj_matrix = symmetric_adj_matrix

    def get_edges_for_protein(self, protein):
        """Gets a list of proteins that the query protein is connected to."""

        protein_index = self.get_protein_index(protein)
        edge_indices = self.adj_matrix[protein_index, :].nonzero()[1].tolist()
        edge_names = [self.get_protein_name_by_index(i) for i in edge_indices]
        return edge_names

    def save_as_pickle(self, save_path):
        """Pickles self to a given file path.

        Parameters
        ----------
        save_path : str
            Desired file path of the output pickle file.
        """
        pickle.dump(self, open(save_path, "wb"))
