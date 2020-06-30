"""
Helper functions for the kinase network diffusion program.

Typical usage example:

    validator = InputValidator()
    legal_ids = validator.validate(list_of_protein_ids)

"""

import pandas as pd

class InputValidator:
    """
    Helper functions for the kinase diffusion app.

    Attributes
    ----------

    uniprot : list
        list of kinase ids in uniprot format (uniprot.com)
    hugo : list
        list of kinase ids in HUGO format (https://www.genenames.org/)
    kinase_fp : str
        file path of the CSV file with the ids
    """

    def __init__(self):
        """
        Inits InputValidator class and loads list of legal kinase ids.
        """
        self.kinase_fp = 'data/list_of_human_kinases.csv'
        self.parse_kinase_file()

    def parse_kinase_file(self):
        """
        Parses uniprot & hugo IDs from the CSV file.
        """
        kinase_df = pd.read_csv(self.kinase_fp)
        self.uniprot = [p.upper() for p in list(kinase_df['uniprot'])]
        self.hugo = [p.upper() for p in list(kinase_df['gene_symbol'])]

    def validate(self, inputs):
        """
        Validates that user supplies protein ids in eitheir uniprot or HUGO format.

        Parameters
        ---------
        imputs : list
            list of protein id strings

        Returns
        -------
        validated : list
            subset of input proteins that are actually present in network
        """

        if not isinstance(inputs, list):
            raise ValueError("expected a list")

        inputs = [protein.upper() for protein in inputs]
        validated = [protein for protein in inputs if protein in self.hugo]
        return validated
