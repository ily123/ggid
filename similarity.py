"""
Calculates GO term similarity between a set of proteins to
produce a similarity matrix, which can then be converted into a
protein-protein network for diffusion.
"""

import multiprocessing
import pickle
import time

import numpy as np
import pandas as pd
from scipy import sparse


class SimilarityCalculator:
    """Calculate Resnik similarity between proteins/entities in a set."""

    def __init__(self, annotations, ontology, ft, proteins, namespace=None):

        self.annotations = annotations.annotations
        self.ontology = ontology
        self.proteins = proteins
        self.ft = ft

        self.mica = {}
        self.ic_mica = {}
        self.counter = 0

        self.proteins = list(set(self.annotations.DB_Object_Symbol) & set(proteins))
        self.proteins.sort()

        if namespace:
            self.annotations = self.annotations.loc[
                self.annotations.Aspect == namespace, :
            ]

        self.annotations2 = self.get_anno_dict(annotations.annotations)

    def sanitize_protein_list(self):
        """Removes proteins with low number of annotations."""

        self.proteins = [p for p in self.proteins if len(self.annotations2[p]) > 9]

    def get_anno_dict(self, annotations):
        """Convert pandas df of entity->term into a dictionary."""

        return dict(annotations.groupby("DB_Object_Symbol")["GO_ID"].apply(list))

    def calculate_similarity(self):
        """Calculate similarity for all-vs-all in the protien set."""

        index_a = []
        index_b = []
        values = []
        for a, protein_a in enumerate(self.proteins):
            for b, protein_b in enumerate(self.proteins):
                if b >= a:
                    sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
                    index_a.append(a)
                    index_b.append(b)
                    values.append(sim)
        sim_matrix = sparse.coo_matrix((values, (index_a, index_b)))
        self.sim_matrix = SimilarityMatrix(sim_matrix, self.proteins)

    def calc_sim_segment(self, list_a, list_all):
        print("worker started")
        t = time.time()
        for protein_a in list_a:
            for protein_b in list_all:
                sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
        print("worker done", time.time() - t)

    def calculate_similarity_mult_cpu(self):
        self.workers = []
        for i in list(range(2)):
            sublist = self.proteins
            p = multiprocessing.Process(
                target=self.calc_sim_segment,
                args=(
                    sublist,
                    self.proteins,
                ),
            )
            self.workers.append(p)
            p.start()

    def calculate_similarity_two_proteins(self, protein_a, protein_b):
        """Calculate Resnik similarity between two proteins."""

        # get protein a terms
        terms_a = self.annotations2[protein_a]
        terms_b = self.annotations2[protein_b]

        # run go terms through a all-vs-all get_mica()
        avg = 0
        best_mica = 0

        row_micas = []
        for term_a in terms_a:
            row_mica = 0
            for term_b in terms_b:
                self.counter += 1
                ic_mica = 0
                if (term_a, term_b) in self.ic_mica:
                    ic_mica = self.ic_mica[(term_a, term_b)]
                else:
                    ic_mica = self.get_ic_mica(term_a, term_b)
                    self.ic_mica[(term_a, term_b)] = ic_mica
                    self.ic_mica[(term_b, term_a)] = ic_mica

                if ic_mica > row_mica:
                    row_mica = ic_mica

                if ic_mica > best_mica:
                    best_mica = ic_mica
            row_micas.append(row_mica)
        #        return best_mica
        #        return avg/(len(terms_a) * len(terms_b))
        return np.array(row_micas).mean()
        #  a b c
        # d x x x
        # e x x x
        # f x x x

    def get_ic_mica(self, term1, term2):
        """Given two terms, Find IC of their  MICA.

        Note:
        IC - information content (frequency of occurance)
        MICA - most informative common ancestor
        """

        ancestors1 = self.ontology.full_ancestry[term1]
        ancestors2 = self.ontology.full_ancestry[term2]

        shared = list(set(ancestors1) & set(ancestors2))
        if len(shared) == 0:
            return 0
        most_specific = self.ft.ic[shared[0]]
        # most specific has smallest absolute IC value
        for ancestor in shared:
            if self.ft.ic[ancestor] > most_specific:
                most_specific = self.ft.ic[ancestor]

        return most_specific

    def get_exhaustive(self):
        """
        Note:
        Going through each protein-protein pair (and the term-term
        pairs within the protein-protien pair) means we will be
        recalculating term-term pairs ad naseum.

        There are n proteins, each has on average 15 terms.

        (15*15 terms) * (19,000 * 19,000) proteins = 90 billion MICA calculations

        If we precalculate the term-term similarity, thats
        47,000 * 47,000 = 2.3 billion MICA calculations
        20,000 * 20,000 proteins = 0.4 billion
        """

        return 0


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

    def enforce_adj_matrix_symmetry():
        """Updates the adjacency matrix to be symmetric about the diagonal."""

        symmetric_adj_matrix = self.adj_matrix + self.adj_matrix.T
        symmetric_adj_matrix[symmetric_adj_matrix > 0] = 1
        self.adj_matrix = symmetric_adj_matrix

    def save_as_pickle(self, save_path):
        """Pickles self to a given file path.

        Parameters
        ----------
        save_path : str
            Desired file path of the output pickle file.
        """

        pickle.dump(self, open(save_path, "wb"))
