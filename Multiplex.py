import pandas as pd



def Check_Multiplex_Compatibility(primer_pairs: pd.DataFrame, heterodimer_max: float = 50.0): # -> pd.DataFrame:
    """
        Check primer sets for multiplex compatibility.
        TODO: Enhance for multiple SNPs.
        - Extend to check cross-SNP interactions (current checks within SNP).
        - Optimize for large primer sets using batch dimer calculations.
        """