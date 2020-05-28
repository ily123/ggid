import time
import multiprocessing
from scipy import sparse

class Similarity:
    """Calculate Resnik similarity between set of proteins/entities"""

    def __init__(self, annotations, ontology, ft, proteins, namespace=None, verbose=False):
        
        self.annotations = annotations.annotations
        self.ontology = ontology
        self.proteins = proteins
        self.ft = ft
        
        self.mica = {}
        self.ic_mica = {}
        self.counter = 0

        #self.proteins = list(set(self.annotations.DB_Object_Symbol) & set(proteins))
        self.proteins = list(set(self.annotations.DB_Object_ID) & set(proteins))
        print(len(self.proteins))
        
        if namespace:
            self.annotations = self.annotations.loc[self.annotations.Aspect == namespace, :]
        
        self.annotations2 = self.get_anno_dict(annotations.annotations)
       
    def get_anno_dict(self, annotations):
        """
        Convert pandas df of entity->term into a dictionary
        of dict[entity]=[term1, term2, termn
        """
        #return dict(annotations.groupby('DB_Object_Symbol')['GO_ID'].apply(list))
        return dict(annotations.groupby('DB_Object_ID')['GO_ID'].apply(list))

    def calculate_similarity(self):
        """Calculate similarity for all-vs-all in the protien set"""
        index_a = []
        index_b = []
        values = []
        for a, protein_a in enumerate(self.proteins):
            for b, protein_b in enumerate(self.proteins):
                if b > a:
                    sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
                    index_a.append(a)
                    index_b.append(b)
                    values.append(sim)
                #print(f'{protein_a} {protein_b} {sim}')
        self.sim_matrix = sparse.coo_matrix((values, (index_a, index_b)))

    def calc_sim_segment(self, list_a, list_all):
        print("worker started")
        t = time.time()
        for protein_a in list_a:
            for protein_b in list_all:
                sim = self.calculate_similarity_two_proteins(protein_a, protein_b)
        print("worker done", time.time()-t)

    def calculate_similarity_mult_cpu(self):
        self.workers = []
        for i in list(range(2)):
            sublist = self.proteins
            p = multiprocessing.Process(target=self.calc_sim_segment, args=(sublist, self.proteins,))
            self.workers.append(p)
            p.start()


    def calculate_similarity_two_proteins(self, protein_a, protein_b):
        """Calculate Resnik similarity between two proteins"""
        
        #get protein a terms
        terms_a = self.annotations2[protein_a]
        terms_b = self.annotations2[protein_b]
        
        # run go terms through a all-vs-all get_mica()
        avg = 0
        for term_a in terms_a:
            for term_b in terms_b:
                self.counter += 1
                if (term_a, term_b) in self.ic_mica:
                    avg += self.ic_mica[(term_a, term_b)]
                else:
                    ic_mica = self.get_ic_mica(term_a, term_b)
                    self.ic_mica[(term_a, term_b)] = ic_mica
                    self.ic_mica[(term_b, term_a)] = ic_mica
                    avg += ic_mica

        return avg/(len(terms_a) * len(terms_b))

    def get_ic_mica(self, term1, term2):
        """
        Given two GO terms, find information content of their
        most informative common ancestor
        """

        ancestors1 = self.ontology.full_ancestry[term1]#.ancestor_ids
        ancestors2 = self.ontology.full_ancestry[term2]#nodes[term2].ancestor_ids
        
        shared = list(set(ancestors1) & set(ancestors2))
        #shared = [val for val in ancestors1 if val in ancestors2]
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
