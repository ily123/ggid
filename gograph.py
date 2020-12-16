import re

import numpy as np
from scipy import sparse


class OBOParser:
    """Gene Ontology .obo file parser."""

    def __init__(self, obo_fp):
        """Inits with file path to obo file."""
        self.obo_fp = obo_fp

    def parse_ontology(self):
        """Parses the .obo file"""
        obo_file = open(self.obo_fp, "r")
        go_dag = GoGraph()
        term_flag = False
        for line in obo_file:
            capture = re.search("(.+):\s(.+)\n", line)
            if line == "[Term]\n":
                term_flag = True
            elif line == "\n":
                term_flag = False

            if capture and term_flag:
                name = capture.group(1)
                content = capture.group(2)
                if name == "id":
                    go_term = GoTerm(id_=content)
                    go_dag.add_node(go_term)  # push term into tree
                elif name in ["is_a", "relationship"]:
                    go_term.is_a.append(content)
                    go_term.ancestor_ids.append(self._parse_ancestor_id(content))
                elif name == "name":
                    go_term.name = content
                elif name == "namespace":
                    go_term.namespace = content
                elif name == "def":
                    go_term.definition = content
        obo_file.close()

        go_dag.draw_connections()  # connect terms
        return go_dag

    @staticmethod
    def _parse_ancestor_id(line):
        """Extracts go term ids from 'is_a' and 'relatioship' fields."""
        capture = re.search("(GO:\d{7})\s!\s", line)
        if capture:
            return capture.group(1)


class GoGraph:
    """Graph implementation of the Gene term Ontology."""

    def __init__(self):
        """Inits naked."""
        # keep track of all nodes using a dictionary
        self.nodes = {}

    def add_node(self, node):
        """Pushes a new node into the node dictionary."""
        self.nodes[node.id] = node

    def get_node_by_id(self, term_id):
        """Returns term node given term id (ex: 'GO:0000001')."""
        if term_id in self.nodes:
            return self.nodes[term_id]

    def draw_connections(self):
        """Connects all ancestor and child nodes."""
        for node in self.nodes.values():
            node_id = node.id
            for ancestor_id in node.ancestor_ids:
                self.connect_two_nodes(node_id, ancestor_id)

    def connect_two_nodes(self, child_id, ancestor_id):
        """Connect a child and an ancestor node given their ids."""
        self.nodes[ancestor_id].add_child(self.nodes[child_id])
        self.nodes[child_id].add_ancestor(self.nodes[ancestor_id])

    def traverse(self, node_id):
        """Traverses GO graph from given node to root.

        Parameters:
            node_id : str
                GO term id of node (ex: 'GO:0000001')
        Return:
            ancestors : list
                list of acestor GO term ids

        Note:
            The list of ancestor ids will have duplicates, because
            in the GO ontology children can have multiple parents and
            so there are multiple paths leading to the root.

            For example:

                    /C\
                A--B   F--root
                    \D/

            Calling traverse(A) will return
                [B, C, F, root, D, F, root]
        """
        ancestors = []
        if self.nodes[node_id].ancestors != []:
            for ancestor in self.nodes[node_id].ancestors:
                ancestors.append(ancestor.id)
                ancestors.extend(self.traverse(ancestor.id))
        return ancestors

    def get_full_ancestry(self, node_id):
        """Returns all ancestors for a given term.

        Parameters:
            node_id : str
                GO term id of node (ex: 'GO:0000001')
        Return:
            ancestors : list
                list of acestor GO term ids
        """
        return list(set(self.traverse(node_id)))

    def make_ancestry_matrix(self):
        """Creates matrix mapping GO terms to their ancestors.

        The ancestry is complete, and includes all terms to the root.
        """
        self.full_ancestry = {}
        for node_id in list(self.nodes):
            self.full_ancestry[node_id] = self.get_full_ancestry(node_id)
        self.ancestry_matrix = AncestryMatrix(self.full_ancestry)


class AncestryMatrix:
    """Helper class to convert ancestry dictionary into a sparse matrix."""

    def __init__(self, full_ancestry):
        """Inits with ancestry dictionary."""
        self._assign_index(full_ancestry)
        self.full_ancestry = full_ancestry
        self.get_matrix()

    def _assign_index(self, full_ancestry):
        """Assign index 0 to N to each term in the dictionary."""
        # do I really need a method for this?
        go_terms = list(full_ancestry)
        indeces = list(range(0, len(go_terms)))
        index_dict = dict(zip(go_terms, indeces))
        self.index_dict = index_dict

    def get_matrix(self):
        """Returns sparse matrix representation of the dictionary."""
        index_a = []
        index_b = []
        for term, ancestors in self.full_ancestry.items():
            term_index = self.index_dict[term]
            for ancestor_term in ancestors:
                index_a.append(term_index)
                index_b.append(self.index_dict[ancestor_term])
        self.matrix = sparse.coo_matrix((np.ones(len(index_a)), (index_a, index_b)))

    def _encode_annotation_counts(self, term_counts):
        """Converts dict of terms counts to a vector.

        If term i is in term_counts, then vector[i] = c,
        and if i is not in term_counts, then vector[i] = 0.
        """
        vector = []
        for term in list(self.index_dict):
            if term in term_counts:
                vector.append(term_counts[term])
            else:
                vector.append(0)
        return np.array(vector).reshape(-1, 1)

    def get_deep_count(self, term_counts):
        """Count implicit number of times a term appears in the ontology."""
        # This part is a little tricky:
        # Whenever a low-level term is assigned to a protein  it's implied that
        # the term's ancestors are assigned to the protein as well. The ancestor
        # terms do not appear in the annotation flat file. So if we want to get
        # the actual number of times each term appears in the corpus, we have to
        # count the number of times it appears implicitly as an ancestror.
        #
        # That is what we do here.
        #
        # We count the number of times each term appears explicitly in the anno
        # coprus, then multiply it by the flat ancestry matrix to crate a matrix
        # of full implied ancestry counts.
        explicit_counts = self._encode_annotation_counts(term_counts)
        return self.matrix.multiply(explicit_counts)


class GoTerm:
    """Node in the Gene Ontology tree/dag."""

    def __init__(self, id_=None, namespace=None, definition=None, is_a=None):
        """Inits empty with a few attrs set to None."""
        self.id = id_
        self.namespace = namespace
        self.definition = definition
        if is_a is None:
            self.is_a = []
        self.ancestor_ids = []
        self.child_ids = []
        self.ancestors = []
        self.children = []

    def add_child(self, child_term):
        """Adds child term."""
        # Why do I have both pointer-like and dictionary tracking?
        self.children.append(child_term)
        self.child_ids.append(child_term.id)

    def add_ancestor(self, ancestor_term):
        """Adds ancestor term."""
        self.ancestors.append(ancestor_term)
