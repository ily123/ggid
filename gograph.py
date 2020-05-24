import numpy as np
import re
from scipy import sparse

class OBOParser:
    """
    Gene Ontology .obo file parser
    """

    def __init__(self, obo_fp):
        self.obo_fp = obo_fp

    def parse_ontology(self):
        """
        Parse the .obo file
        """

        obo_file = open(self.obo_fp, 'r')
        go_list = {}
        blank_template = {'id': None, 'name': None, 'namespace': None,
                          'def': None, 'is_a': []}

        go_dag = GoGraph()
        term_flag = False
        for line in obo_file:
            capture = re.search('(.+):\s(.+)\n', line)
            if line == '[Term]\n':
                term_flag = True
            elif line == '\n':
                term_flag = False
                
            if capture and term_flag:
                name = capture.group(1)
                content = capture.group(2)
                if name == 'id':
                    go_term = GoTerm(id_=content)
                    go_dag.add_node(go_term) # push term into tree
                elif name == 'is_a' or name == 'relationship':
                    go_term.is_a.append(content)
                    go_term.ancestor_ids.append(self.parse_ancestor_id(content))
                elif name == 'name':
                    go_term.name = content
                elif name == 'namespace':
                    go_term.namespace = content
                elif name == 'def':
                    go_term.definition = content
        obo_file.close()

        go_dag.draw_connections() # connect terms
        return go_dag

    def parse_ancestor_id(self, line):
        """
        Extracts go term ids from 'is_a' and 'relatioship' fields
        """
        capture = re.search('(GO:\d{7})\s!\s', line)
        if capture:
            return capture.group(1)



class GoGraph:
    """
    Graph implementation of the Gene term Ontology

    Note:

        Implementation closely follows the one described in chapter 8
        of Runestone's academy data strcutres course

        https://runestone.academy/runestone/books/published/
            pythonds/Graphs/Implementation.html
    """


    def __init__(self):
        """
        Keep track of all nodes using a dictionary
        """
        
        self.nodes = {}

    def add_node(self, node):
        """
        Push a new node into the node dictionary
        """
        
        self.nodes[node.id] = node

    def get_node(self, term_id):
        """
        Return term node given term id (ex: 'GO:0000001')
        """

        if term_id in self.nodes:
            return self.nodes[term_id]
        else:
            return None
    
    def draw_connections(self):
        """
        Connect all ancestor and child nodes
        """

        for node in self.nodes.values():
            node_id = node.id
            for ancestor_id in node.ancestor_ids:
                self.connect_two_nodes(node_id, ancestor_id)

    def connect_two_nodes(self, child_id, ancestor_id):
        """
        Connect a child and an ancestor node given their term ids

        """

        self.nodes[ancestor_id].add_child(self.nodes[child_id])
        self.nodes[child_id].add_ancestor(self.nodes[ancestor_id])
         
    def traverse(self, node_id):
        """
        Traverse GO graph from given node to root
        
        Attributes:
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
        """
        Get all ancestors for a given term
        
        Attributes:
            node_id : str
                GO term id of node (ex: 'GO:0000001')
        Return:
            ancestors : list
                list of acestor GO term ids
        """

        return list(set(self.traverse(node_id)))

    def make_ancestry_matrix(self):
        """
        Create matrix mapping GO terms to their complete chain
        of ancestor terms, up to the root
        """
        self.full_ancestry = {}
        for node_id in list(self.nodes):
            self.full_ancestry[node_id] = self.get_full_ancestry(node_id)
    
        
        self.ancestry_matrix = AncestryMatrix(self.full_ancestry)
        

class AncestryMatrix:
    """
    Helper class to convert ancestry dictionary into a sparse matrix

    Note:
        This class is meant to simplify my life, as I don't want to write
        a bunch of for-loops. It's probably a speed up as well.
    """

    def __init__(self, full_ancestry):
        self.assign_index(full_ancestry)
        self.full_ancestry = full_ancestry
        self.get_matrix()

    def assign_index(self, full_ancestry):
        go_terms = list(full_ancestry)
        indeces = list(range(0,len(go_terms)))
        index_dict = dict(zip(go_terms, indeces))
        self.index_dict = index_dict

    def get_matrix(self):
        index_a = []
        index_b = []

        for term, ancestors in self.full_ancestry.items():
            term_index = self.index_dict[term]
            for ancestor_term in ancestors:
                index_a.append(term_index)
                index_b.append(self.index_dict[ancestor_term])
        
        self.matrix = sparse.coo_matrix((np.ones(len(index_a)), (index_a, index_b))) 

    def encode_annotation_vector(self, term_counts):
        
        vector = [] 
        for term in list(self.index_dict):
            if term in term_counts:
                vector.append(term_counts[term])
            else:
                vector.append(0)
        
        return np.array(vector).reshape(-1,1)       

    def get_deep_count(self, term_counts):

        vector = self.encode_annotation_vector(term_counts)
        return self.matrix.multiply(vector)


class GoTerm:
    """
    Node in the Gene Ontology tree/dag
   
    """

    def __init__(self, id_=None, namespace=None, definition=None, is_a=None):

        self.id = id_
        self.namespace = namespace
        self.definition = definition
        if is_a == None:
            self.is_a = []

        self.ancestor_ids = []
        self.child_ids = []
        self.ancestors = []
        self.children = []
        
    def add_child(self, child_term):
        self.children.append(child_term)
        self.child_ids.append(child_term.id)

    def add_ancestor(self, ancestor_term):
        self.ancestors.append(ancestor_term)

