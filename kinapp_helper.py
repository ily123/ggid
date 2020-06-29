import pandas as pd

class InputValidator:
    """
    Helper functions for the kinase diffusion app

    Attributes
    ----------

    uniprot : list
        list of kinase ids in uniprot format (uniprot.com)
    hugo : list
        list of kinase ids in HUGO format (https://www.genenames.org/)
    fp : str
        file path of the CSV file with the ids
    """

    def __init__(self):

        self.fp = 'data/list_of_human_kinases.csv'
        self.parse_kinase_file()
    
    def parse_kinase_file(self):
        """
        Parses uniprot & hugo IDs from the CSV file
        """
        kinase_df = pd.read_csv(self.fp)
        self.uniprot = [p.upper() for p in list(kinase_df['uniprot'])]
        self.hugo = [p.upper() for p in list(kinase_df['gene_symbol'])]
    
    def validate(self, inputs):
        """
        Validates that user supplies protein ids in eitheir uniprot or HUGO format

        Parameters
        ---------
        imputs : list
            list of protein id strings

        Returns
        -------
        validated : list
            list of input proteins that are present in network
        """

        assert type(inputs) == list, "expected a list"
        
        inputs = [p.upper() for p in inputs]
        validated = []
        for protein in inputs:
            if protein in self.hugo:
                validated.append(protein)
        return validated
    
    def blah(self):
        print(0)
