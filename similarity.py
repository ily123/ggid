import time

class Similarity:
    """Calculate Resnik similarity between set of proteins/entities"""

    def __init__(self, annotations, ontology, proteins, verbal=False):
        
       # if verbal:
       #     t0 = time.time()
       #     print('Parsing ontology file')
       #     self.annotations = annotations
       #     print('....done in %1.1f sec' % (time.time() - t0))
       #     t0 = time.time()
       #     self.ontology = ontology
       #     self.proteins = proteins
       # else:
       self.annotations = annotations
       self.ontology = ontology
       self.proteins = proteins

    def calculate_similarity(self):
        """Calculate similarity for all-vs-all in the protien set"""
        return 0

    def calculate_similarity_two_proteins(self, protein_a, protein_b):
        """Calculate Resnik similarity between two proteins"""
        
        #get protein a terms

        # get protein b terms
        
        # run go terms through a all-vs-all get_mica()
        return 0

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

    def get_ic_mica(self, term1, term2):
        """
        Given two GO terms, find information content of their
        most informative common ancestor
        """

        ancestors1 = term1.ancestors
        ancestors2 = term2.ancestors
        
        shared = list(set(ancestors1) & set(ancestors2))
        most_specific = shared[0]
        # most specific has smallest absolute IC value
        for ancestor in shared:
            if ancestor.get_ic < most_specific.get_ic:
                most_specific = ancestor.get_ic
        
        return 0

