import time

class Similarity:
    """Calculate Resnik similarity between set of proteins/entities"""

    def __init__(self, annotations, ontology, ft, proteins, namespace=None, verbose=False):
        
        self.annotations = annotations.annotations
        self.ontology = ontology
        self.proteins = proteins
        self.ft = ft
        
        self.mica = {}
        
        if namespace:
            self.annotations = self.annotations.loc[self.annotations.Aspect == namespace, :]
        
    def calculate_similarity(self):
        """Calculate similarity for all-vs-all in the protien set"""
        for protein_a in self.proteins:
            for protein_b in self.proteins:
                sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
                #print(f'{protein_a} {protein_b} {sim}')

    def calculate_similarity_two_proteins(self, protein_a, protein_b):
        """Calculate Resnik similarity between two proteins"""
        
        #get protein a terms
        terms_a = list(self.annotations.loc[self.annotations.DB_Object_Symbol==protein_a, 'GO_ID'])
        terms_b = list(self.annotations.loc[self.annotations.DB_Object_Symbol==protein_b, 'GO_ID'])
        
        # run go terms through a all-vs-all get_mica()
        avg = 0
        for term_a in terms_a:
            for term_b in terms_b:
                avg += self.get_ic_mica(term_a, term_b)
        return avg/(len(terms_a) * len(terms_b))

    def get_ic_mica(self, term1, term2):
        """
        Given two GO terms, find information content of their
        most informative common ancestor
        """

        ancestors1 = self.ontology.full_ancestry[term1]#.ancestor_ids
        ancestors2 = self.ontology.full_ancestry[term2]#nodes[term2].ancestor_ids
        
        shared = list(set(ancestors1) & set(ancestors2))
        if len(shared) == 0:
            return 0
       # print('---')
       # print(term1)
       # print(ancestors1)
       # print(term2)
       # print(ancestors2)
       # print(shared)
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
