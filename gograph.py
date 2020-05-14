
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
        Return term node given term id (ex: GO:0000001)
        """

        if term_id in self.nodes:
            return self.nodes[term_id]
        else:
            return None

    def connect_nodes(self, ancestor_id, child_id):
        """
        Connect a child and an ancestor node given term ids

        Note: nodes must already be present in self.nodes dict
        """
        
        self.nodes[ancestor_id].set_child(self.nodes[child_id])
        self.nodes[child_id].set_ancestor(self.nodes[ancestor_id])
        

class GoNode:
    """
    fff
    """
    
    def __init__(self, payload)
        self.payload = payload
        self.ancestors = []
        self.children = []

    def set_child(self, node):
        self.child.append(node)
    
    def set_ancestor(self, node):
        self.ancestors.append(node)


class GoTerm:
    """
    Payload that goes with GoNode
   
    This contains all the pertinent information about the term
    """


   # insert parsing logic here 
