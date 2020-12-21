import re
from typing import Any, Dict, List, Type, Union

import numpy as np
import pandas as pd
from scipy import sparse


class GoTerm:
    """GO term node in the Gene Ontology tree/dag.

    Attributes
    ----------
    id_ : str
        GO id of the term
    namespace : str, optional
        ontology namespace of the term (MF, BP, or CC)
    definition : str, optional
        definition of the term
    is_a : str, optional
        is_a ancestor field
    ancestors : List[GoTerm]
        list of the term's ancestor terms
    children : list[GoTerm]
        list of the term's child terms
    """

    def __init__(self, id_: str) -> None:
        """Inits with GO term id.

        Parameters
        ----------
        id_ : str
            GO term id (ex: "GO:0000001")
        """
        self.id = id_
        self.namespace = None
        self.definition = None
        self.is_a = None
        self.ancestors = []
        self.children = []

    def add_child(self, child_term: Type["GoTerm"]) -> None:
        """Adds child term.

        Parameters
        ----------
        child_term : GoTerm
            GoTerm object that is a child of this term
        """
        self.children.append(child_term)

    def add_ancestor(self, ancestor_term: Type["GoTerm"]) -> None:
        """Adds ancestor term.

        Parameters
        ----------
        ancestor_term : GoTerm
            GoTerm object that is an ancestor of this term
        """
        self.ancestors.append(ancestor_term)


class GoGraph:
    """Graph implementation of the Gene term Ontology.

    Attributes
    ----------
    nodes : Dict
        dictionary of nodes in the graph, addressed by id
    full_ancestry : Dict[List]
        dictionary of term -> term ancestors traced to the root
    """

    def __init__(self) -> None:
        """Inits empty."""
        # keep track of all nodes using a dictionary
        self.nodes = {}

    def add_node(self, node: Type[GoTerm]) -> None:
        """Pushes a new node into the node dictionary.

        Parameters
        ----------
        node : GoTerm
            a GO term object
        """
        self.nodes[node.id] = node

    def get_node_by_id(self, term_id: str) -> Union[None, Type[GoTerm]]:
        """Returns term node given term id (ex: 'GO:0000001').

        Parameters
        ----------
        term_id : str
            GO id of the term (ex: "GO:0000001")
        """
        if term_id in self.nodes:
            return self.nodes[term_id]
        return None

    def draw_connections(self) -> None:
        """Connects all ancestor and child nodes."""
        for node in self.nodes.values():
            node_id = node.id
            for ancestor_id in node.ancestor_ids:
                self.connect_two_nodes(node_id, ancestor_id)

    def connect_two_nodes(self, child_id: str, ancestor_id: str) -> None:
        """Connects a child and an ancestor node given their ids.

        Parameters
        ----------
        child_id : str
            GO id of the child term
        ancestor_id : str
            GO id of the ancestor term
        """
        self.nodes[ancestor_id].add_child(self.nodes[child_id])
        self.nodes[child_id].add_ancestor(self.nodes[ancestor_id])

    def traverse(self, node_id: str) -> List[str]:
        """Traverses GO graph from given node to root.

        Parameters
        ----------
        node_id : str
            GO term id of node (ex: 'GO:0000001')

        Returns
        -------
        ancestors : list
            list of acestor GO term ids

        Note:
            The list of ancestor ids will have duplicates, because
            in the GO ontology children can have multiple parents and
            so there are multiple paths leading to the root.

            For example:

                    -C-
                A--B   F--root
                    -D-

            Calling traverse(A) will return
                [B, C, F, root, D, F, root]
        """
        ancestors = []
        if self.nodes[node_id].ancestors != []:
            for ancestor in self.nodes[node_id].ancestors:
                ancestors.append(ancestor.id)
                ancestors.extend(self.traverse(ancestor.id))
        return ancestors

    def get_full_ancestry(self, node_id: str) -> List[str]:
        """Returns all ancestors for a given term.

        Parameters
        ----------
            node_id : str
                GO term id of node (ex: 'GO:0000001')
        Returns
        -------
            ancestors : list
                list of acestor GO term ids
        """
        return list(set(self.traverse(node_id)))

    def make_ancestry_matrix(self) -> None:
        """Creates matrix mapping GO terms to their ancestors.

        The ancestry is complete, and includes all terms to the root.
        """
        self.full_ancestry = {}
        for node_id in list(self.nodes):
            self.full_ancestry[node_id] = self.get_full_ancestry(node_id)
        self.ancestry_matrix = AncestryMatrix(self.full_ancestry)


class Annotations:
    """Container for storing annotation data."""

    def __init__(self, anno_fp: str) -> None:
        """Inits with path to annotations gaf file.

        Parameters
        ----------
        anno_fp : str
            path to annotations file in gaf-2 format
        """
        self.anno_fp = anno_fp
        header = [
            "DB",
            "DB_Object_ID",
            "DB_Object_Symbol",
            "Qualifier",
            "GO_ID",
            "DB_Reference",
            "Evidence_Code",
            "With_or_From",
            "Aspect",
            "DB_Object_Name",
            "DB_Object_Synonym",
            "DB_Object_Type",
            "Taxon",
            "Date",
            "Assigned_By",
            "Annotation_Extension",
            "Gene_Product_Form_ID",
        ]
        self.annotations = pd.read_csv(anno_fp, header=None, sep="\t", comment="!")
        self.annotations.columns = header

    def get_counts(self) -> Dict[str, int]:
        """Counts number of times terms appear in annotation corpus.

        Returns
        -------
        counts : Dict[str, int]
            term counts
        """
        counts = self.annotations.groupby("GO_ID").size()
        return counts

    def filter(self, keep_namespace: str) -> None:
        """Removes annotations outside of namespace, inplace.

        Parameters
        ----------
        keep_namespace : str
            string token of namespace to keep (C, P, or F)
        """
        legal_namespace = ["C", "P", "F"]
        if keep_namespace.upper() not in legal_namespace:
            raise ValueError(
                "Namespace must be one of: %s" % ", ".join(legal_namespace)
            )
        self.annotations = self.annotations.loc[
            self.annotations.Aspect == keep_namespace, :
        ]


class FrequencyTable:
    def __init__(self, annotations, ancestry_matrix):
        shallow_count = annotations.get_counts()
        self.deep_count = ancestry_matrix.get_deep_count(shallow_count)
        self.index_dict = ancestry_matrix.index_dict

        self.information_content = -1 * np.log10(
            self.deep_count.sum(axis=0) / self.deep_count.sum()
        )
        self.ic = dict(zip(list(self.index_dict), self.information_content.tolist()[0]))


class OBOParser:
    """Gene Ontology .obo file parser."""

    def __init__(self, obo_fp: str) -> None:
        """Inits with file path to obo file."""
        self.obo_fp = obo_fp

    def parse_ontology(self) -> Type[GoGraph]:
        """Parses the .obo file."""
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
    def _parse_ancestor_id(line) -> Union[None, str]:
        """Extracts go term ids from 'is_a' and 'relatioship' fields."""
        capture = re.search("(GO:\d{7})\s!\s", line)
        if capture:
            return capture.group(1)


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
        # coprus, then multiply it by the flat ancestry matrix to create a matrix
        # of complete/implied ancestry counts.
        explicit_counts = self._encode_annotation_counts(term_counts)
        return self.matrix.multiply(explicit_counts)
