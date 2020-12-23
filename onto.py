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
        self.name = None
        self.definition = None
        self.ancestor_ids = []
        self.ancestors = []
        self.children = []
        self.specificity = None

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
    """

    def __init__(self, obo_fp: str) -> None:
        """Inits with obo file fp.

        Parameters
        ----------
        obo_fp : str
            file path to the OBO file containing the GO DAG
        """
        # keep track of all nodes using a dictionary
        self.nodes = {}
        self.obo_fp = obo_fp

    def parse_ontology(self) -> None:
        """Parses the .obo file."""
        obo_file = open(self.obo_fp, "r")
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
                    self.add_node(go_term)  # push term into dict
                elif name in ["is_a", "relationship"]:
                    capture = re.search("(GO:\d{7})\s!\s", content)
                    ancestor_id = capture.group(1)
                    go_term.ancestor_ids.append(ancestor_id)
                elif name == "name":
                    go_term.name = content
                elif name == "namespace":
                    go_term.namespace = content
                elif name == "def":
                    go_term.definition = content
                elif name == "is_obsolete":
                    # obsolete status is the last line in the term entry
                    # we don't want obsolete terms in the tree
                    del self.nodes[go_term.id]
        obo_file.close()
        self.draw_connections()  # connect terms

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

    def calculate_term_specificity(
        self, annotations: Type["Annotations"]
    ) -> Dict[str, float]:
        """Calculates specificity of each GO term in the graph.

        Parameters
        ----------
        annotations : Annotations
            GO term annotation corpus, a mapping of annotations to terms

        Returns
        -------
        term_specificity : dict[str, float]
            mapping of terms to their specificty (neg log frequency)
        """
        full_ancestry = {}
        for node_id in list(self.nodes):
            full_ancestry[node_id] = self.get_full_ancestry(node_id)
        spec_calc = SpecificityCalculator(full_ancestry, annotations)
        term_specificity = spec_calc.get_specificity()
        return term_specificity

    def assign_term_specificity(self, term_specificity: Dict[str, float]) -> None:
        """Assigns term specificity to each term in the tree.

        Parameters
        ----------
        term_specificity : dict[str, float]
            mapping of terms to their specificty (neg log frequency)
        """
        for term, specificity in term_specificity.items():
            self.nodes[term].specificity = specificity


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

    def get_as_dict(self) -> Dict[str, List[str]]:
        """Returns annotations as protein->terms dict."""
        return dict(self.annotations.groupby("DB_Object_Symbol")["GO_ID"].apply(list))


class SpecificityCalculator:
    """Helper class to calculate frequency (specificity) of terms."""

    def __init__(
        self, full_ancestry: Dict[str, List[str]], annotation_corpus: Type[Annotations]
    ) -> None:
        """Inits with ancestry dictionary.

        Parameters
        ----------
        full_ancestry : Dict[str, List[str]]
            Dictionary mapping term ids to ancestor ids (traced to root)
        annotation_corpus: onto.Annotations object
            Dictionary mapping terms to their counts in the anno corpus
        """
        self.full_ancestry = full_ancestry
        self.annotation_corpus = annotation_corpus
        self._ancestry_matrix = self._get_matrix()

    def _get_matrix(self) -> Type[sparse.coo_matrix]:
        """Returns sparse matrix representation of the ancestry dictionary.

        Returns
        -------
        ancestry_matrix : sparse.coo_matrix
            matrix where cell i, j is set to 1 if terms i, j are a child-ancestor pair
        """
        index_a = []
        index_b = []
        term_index = dict(
            zip(sorted(list(self.full_ancestry)), range(len(self.full_ancestry)))
        )
        for term, ancestors in self.full_ancestry.items():
            for ancestor_term in ancestors:
                index_a.append(term_index[term])
                index_b.append(term_index[ancestor_term])
        ancestry_matrix = sparse.coo_matrix((np.ones(len(index_a)), (index_a, index_b)))
        return ancestry_matrix

    def _encode_annotation_counts(self) -> np.array:
        """Creates a vector of shallow term counts.

        For vector[i] = c, c is the number of times term i appears in the anno corpus.
        """
        term_counts = self.annotation_corpus.get_counts()
        vector = []
        for term in sorted(list(self.full_ancestry)):
            if term in term_counts:
                vector.append(term_counts[term])
            else:
                vector.append(0)
        return np.array(vector).reshape(-1, 1)

    def get_full_count(self) -> np.array:
        """Count full number of times a term appears in the ontology.

        Returns
        -------
        full_count : np.array
            array full_count[i] = c, where c is the number of times term i appears
            in the anno corpus, plus all the number of times it's implied as an ancestor
        """
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
        explicit_counts = self._encode_annotation_counts()
        full_count = self._ancestry_matrix.multiply(explicit_counts).sum(axis=0)
        return full_count

    def get_specificity(self) -> Dict[str, float]:
        """Calculates specificity of each term using full counts.

        Returns
        -------
        term_specificity : Dict[str, float]
            dictionary mapping terms to their negative log frequency (specificity)
        """
        full_count = self.get_full_count()
        specificity = -1 * np.log10(full_count / full_count.sum())
        term_specificity = dict(
            zip(sorted(list(self.full_ancestry)), specificity.tolist()[0])
        )
        return term_specificity
