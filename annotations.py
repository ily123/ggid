import numpy as np
import pandas as pd

class Annotations:

    def __init__(self, anno_fp):
        self.anno_fp = anno_fp
        header = [
            'DB',
            'DB_Object_ID',
            'DB_Object_Symbol',
            'Qualifier',
            'GO_ID',
            'DB_Reference',
            'Evidence_Code',
            'With_or_From',
            'Aspect',
            'DB_Object_Name',
            'DB_Object_Synonym',
            'DB_Object_Type',
            'Taxon',
            'Date',
            'Assigned_By',
            'Annotation_Extension',
            'Gene_Product_Form_ID'
            ]
        self.annotations = pd.read_csv(anno_fp, header = None,
                sep= '\t', comment = '!')
        self.annotations.columns = header
        
    def parse(self):
        """
        Parse GO annotation file (.gaf)
        """
        return 0

    def get_counts(self):
        """
        Count number of times each GO term appears
        in annotation corpus
        """
        counts = self.annotations.groupby('GO_ID').size()
        return counts
    
class FrequencyTable():

    def __init__(self, annotations, ancestry_matrix):
        shallow_count = annotations.get_counts()
        self.deep_count = ancestry_matrix.get_deep_count(shallow_count)
        self.index_dict = ancestry_matrix.index_dict

        self.information_content = -1*np.log(self.deep_count.sum(axis=0)/self.deep_count.sum())
        self.ic = dict(zip(list(self.index_dict),
            self.information_content.tolist()[0]))
