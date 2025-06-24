import pandas as pd


# Sequence -> Sequence
def Reverse_Compliment(sequence: str): #  -> str
    # return str(Seq.Seq(sequence).reverse_complement())
    ...


def Introduce_Mismatch(sequence: str, position: int): #  -> List[str]:
    """
        Introduce a mismatch at the antepenultimate position (3rd from last).
        TODO: Validate mismatch rule for biological accuracy.
        - Consider additional mismatch types if needed for specificity.
        - Add checks for invalid sequences (e.g., non-ACGT bases).
        """


def Evaluate_Primers(primer_seq: str): # -> Dict:
    """
        Evaluate primer quality using primer3-py.
        TODO: Enhance error handling.
        - Handle primer3-py failures gracefully.
        - Add logging for failed evaluations.
        """
    ...

def Rank_Primers(primers: pd.DataFrame): # -> pd.DataFrame:
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    ...

def Filter_Primers(primers: pd.DataFrame, tm_min: float = 60.0, tm_max: float = 65.0, hairpin_max: float = 45.0, homodimer_max: float = 45.0): # -> pd.DataFrame:
    """
        Filter primers based on quality metrics.
        TODO: Add fallback for strict filtering.
        - If no primers pass, consider relaxing criteria or selecting best available.
        """


def Generate_Allele_Spesific_Primers(snp_data: pd.DataFrame, min_len: int = 18, max_len: int = 28): # -> pd.DataFrame:
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        """
    # make an empty list of primers
    #iterate through each function calling the find primers function. 
    #this function is a paralizing shell that will house the true find primers function.

def Find_Primers(snp_row, min_len, max_len):
        this_allele_primers = []
        snp_id = snp_row["snpID"]
        allele = snp_row["allele"]
        sequence = snp_row["sequence"]
        center = snp_row["position"]

        forward = sequence[center - max_len + 1:center + 1]#this gets the largest segment.
        forward_mismatch = introduce_mismatch(forward, max_len)
        sequence_length = len(forward_mismatch)

        if sequence_length >= min_len:
            for length in range(max_len-min_len-1):#possibe bug if the forward missmatch is smaller than the minimum length
                trimmed = forward_mismatch[length]
                this_allele_primers.append({
                    "snpID": snp_id,
                    "allele": allele,
                    "primer_sequence": trimmed,
                    "direction": "forward",
                    "length": sequence_length-length
                })
        else:
            print("there's a troll in the dungeon!!!")


            #do the rest later
        # Reverse primer: downstream sequence, reverse complemented
        reverse = reverse_complement(sequence[center:center + length])
        reverse_mismatches = introduce_mismatch(reverse, length)
        for rm in reverse_mismatches:
            primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": rm,
                "direction": "reverse",
                "length": length
            })
    return pd.DataFrame(primers)

def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """