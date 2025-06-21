# from Multiplex import *
# from Output import *
# from Primer_functions import *
# import primer3
import pcr_lib as pl



testvar = "ATGA"

print(pl.reverse_complement(testvar))

test_var = pl.snp_data()
print(test_var)
test_var.allele = "awiohduiahwuidhuiawhuid"
print(test_var.allele)
print(test_var["allele"])



def Fetch_SNP_Data(rsids: list[str], flank_length: int = 800): #  -> pd.DataFrame
    """
        Fetch SNP data from dbSNP (mocked for demo).
        TODO: Replace with real Entrez/API call to fetch rsIDs, alleles, and flanking sequences from dbSNP.
        - Use Bio.Entrez.esearch and efetch to query SNP database.
        - Parse XML/JSON output for alleles and flanking sequences.
        - Handle errors for invalid rsIDs or missing data.
        """




def Main():
    """
        Main function to generate, filter, rank, pair, and export primers.
        TODO: Add comprehensive error handling and logging.
        - Log progress and errors to a file.
        - Add input validation for rsIDs and output format.
        TODO: Test with real SNP data.
        - Validate output with biological experts.
        - Benchmark performance for large SNP sets.
        """

