import pandas as pd
import re 
from Bio.Seq import Seq

# # Sequence -> Sequence
# def Reverse_Complement(sequence: str): #  -> str
    
#     return str(Seq(sequence).reverse_complement())
# This function was unnecessary


def Introduce_Mismatch(primer_sequence: str) -> str:
    """
    Introduces a base mismatch at the antepenultimate position (3rd from last).
    """
    # Ensure valid string input
    if not primer_sequence or not isinstance(primer_sequence, str):
        print("Warning: Invalid primer input.")
        return primer_sequence

    primer_sequence = primer_sequence.upper().strip()

    # Must only contain A, C, G, T
    if not re.match("^[ACGT]+$", primer_sequence):
        print(f"Warning: Invalid characters in primer: {primer_sequence}")
        return primer_sequence

    # Must be long enough to have a 3rd-to-last base
    if len(primer_sequence) < 3:
        print(f"Warning: Primer too short for mismatch: {primer_sequence}")
        return primer_sequence

    # Simple mismatch rules (purine↔purine, pyrimidine↔pyrimidine)
    mismatch_rules = {
        "A": "G", "G": "A",
        "C": "T", "T": "C"
    }

    pos = len(primer_sequence) - 3  # Antepenultimate index
    base = primer_sequence[pos]
    mismatch = mismatch_rules.get(base)

    if mismatch is None:
        print(f"Warning: No mismatch rule for base '{base}'")
        return primer_sequence

    # Replace the base with its mismatch
    return primer_sequence[:pos] + mismatch + primer_sequence[pos + 1:]


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
    all_primers = []
    #made these two separate functions incase we want to parallelize them
    for _, row in snp_data.iterrows():#left this as it was. Need to look at how the data will actually come in
        all_primers.append(Find_Primers(row, min_len, max_len))#call the function again and again. My say is we batch by snpID or something and multiprocess a batch
        #starting a whole thread just for a couple for loops doesn't quite seem justified.

    # make an empty list of primers
    #iterate through each function calling the find primers function. 
    #this function is a paralizing shell that will house the true find primers function.

def Find_Primers(snp_row, min_len, max_len):
    this_allele_primers = []
    snp_id = snp_row["snpID"]
    allele = snp_row["allele"]
    sequence = snp_row["sequence"]
    center = snp_row["position"]
    print(f"This is the sequence: {sequence}")

    forward = sequence[center - max_len + 1 :center + 1]#this gets the largest segment.
    print(f"This is the forward: {forward}")
    forward_mismatch = Introduce_Mismatch(forward)
    print(f"this is the forward after the mismatch: {forward_mismatch}")
    forward_length = len(forward_mismatch)
    print(f" this is it's length: {forward_length}")
    print(f" this is the minimum length: {min_len}")

    if forward_length >= min_len:
        for length in range(max_len-(min_len-1)):#possibe bug if the forward missmatch is smaller than the minimum length
            trimmed = forward_mismatch[length:]
            print(trimmed)
            this_allele_primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": "forward",
                "length": forward_length-length
            })
            
    else:
        print("there's a troll in the dungeon!!!")


    # Reverse primer: downstream sequence, reverse complemented
    # print(f"here's the sequence again: {sequence}")
    reverse = str(Seq(sequence[center + 1:center+max_len-1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
    print(f"the reverse re-stringafied: {reverse}")
    reverse_mismatch = Introduce_Mismatch(reverse)
    print(f"the missmatch reverse: {reverse_mismatch}")
    reverse_length = len(reverse_mismatch)
    print(f" and the mismatch reverse length: {reverse_length}")

    if reverse_length >= min_len:
        for length in range(max_len-(min_len-1)):
            trimmed = forward_mismatch[length:]
            print(trimmed)
            this_allele_primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": "reverse",
                "length": reverse_length-length
            })
            
    else:
        print("there's a troll in the dungeon!!! (reverse)")
    return pd.DataFrame(this_allele_primers)

def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """