def Reverse_Compliment():
    # return str(Seq.Seq(sequence).reverse_complement())
    ...


def Introduce_Mismatch():
    """
        Introduce a mismatch at the antepenultimate position (3rd from last).
        TODO: Validate mismatch rule for biological accuracy.
        - Consider additional mismatch types if needed for specificity.
        - Add checks for invalid sequences (e.g., non-ACGT bases).
        """


def Evaluate_Primers():
    """
        Evaluate primer quality using primer3-py.
        TODO: Enhance error handling.
        - Handle primer3-py failures gracefully.
        - Add logging for failed evaluations.
        """
    ...

def Rank_Primers():
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    ...

def Filter_Primers():
    """
        Filter primers based on quality metrics.
        TODO: Add fallback for strict filtering.
        - If no primers pass, consider relaxing criteria or selecting best available.
        """


def Generate_Allele_Spesific_Primers():
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        """
    ...

def Generate_Matching_Primers():
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """