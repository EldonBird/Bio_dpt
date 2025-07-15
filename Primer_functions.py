import pandas as pd
import re 
from Bio.Seq import Seq
import logging
from typing import Dict
import primer3

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


def Evaluate_Primers(primer_dict: dict): # -> Dict:
    """
        Evaluate primer quality using primer3-py.
        TODO: Enhance error handling.
        - Handle primer3-py failures gracefully.
        - Add logging for failed evaluations.
        """

    # Setup logging
    logging.basicConfig(
        # print('in set up logic')
        filename="primer_evaluation.log",
        level=logging.INFO,
        # print('logging info')
        format="%(asctime)s - %(levelname)s - %(message)s"
    )


    try:

        # Use primer3's analysis functions instead of design_primers
        tm_result = primer3.bindings.calc_tm(primer_dict['primer_sequence'])
        # gc_content = primer3.bindings.calc_gc(primer_seq)
        hairpin = primer3.bindings.calc_hairpin(primer_dict['primer_sequence'])
        homodimer = primer3.bindings.calc_homodimer(primer_dict['primer_sequence'])

        primer_dict["tm"] = tm_result
        # "gc_content": gc_content,
        primer_dict["hairpin"] = hairpin.tm if hasattr(hairpin, 'tm') else 0
        primer_dict["homodimer"] = homodimer.tm if hasattr(homodimer, 'tm') else 0
        primer_dict["success"] = True
        return primer_dict
    except Exception as e:
        logging.error(f"Primer evaluation failed for sequence: {primer_dict['primer_sequence']}\nError: {str(e)}")
        primer_dict["tm"] = 0
        primer_dict["hairpin"] = 999
        primer_dict["homodimer"] = 999
        primer_dict["success"] = False
        return primer_dict

def rank_primers(primers: list[dict], target_tm = 62.5, target_gc = 50) -> pd.DataFrame:
    """
    Rank primers based on Tm proximity to 62.5Â°C and GC content.
    TODO: Refine ranking criteria.
    - Consider weighting Tm vs. GC scores.
    - Add user-configurable ranking metrics.
    """
    for primer in primers:
        primer = Evaluate_Primers(primer)
        primer["tm_score"] = abs(primer["tm"] - target_tm)
        # primer["gc_score"] = abs(primer["gc_content"] - target_gc)
        primer["score"] = primer["tm_score"] + primer["hairpin"] + primer["homodimer"]#+ primer["gc_score"] 

    return primers.sort_values("score").groupby(["snpID", "allele", "direction"]).head(5)

    

def Filter_Primers(primers: pd.DataFrame, tm_min: float = 60.0, tm_max: float = 65.0, hairpin_max: float = 45.0, homodimer_max: float = 45.0): # -> pd.DataFrame:
    """
        Filter primers based on quality metrics.
        TODO: Add fallback for strict filtering.
        - If no primers pass, consider relaxing criteria or selecting best available.
        """


def Generate_Allele_Specific_Primers(snps_list: list[dict], min_len: int = 18, max_len: int = 28): # -> pd.DataFrame:
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        """
    all_primers = []
    min_len -= 2 #don't know why but these need to have 2 and 1 removed from the inputs to get the desired lengths
    max_len -= 1

    for snp_dict in snps_list:
        this_allele_primers = []#a list of dictionaries
        snp_id = snp_dict["snpID"]
        allele = snp_dict["allele"]
        sequence = snp_dict["sequence"]
        snp_pos = snp_dict["position"]

        forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
        forward_mismatch = Introduce_Mismatch(forward)
    
        reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
        reverse_mismatch = Introduce_Mismatch(reverse)

        this_allele_primers.extend(Make_Primers(this_allele_primers, forward_mismatch, min_len, max_len, snp_id, allele))
        this_allele_primers.extend(Make_Primers(this_allele_primers, reverse_mismatch, min_len, max_len, snp_id, allele, "reverse"))
 
        all_primers.append(this_allele_primers)
    return all_primers
   

def Make_Primers(list, seq, min_len, max_len, snp_id, allele, direction="forward"):
    seq_length = len(seq)

    if seq_length >= min_len:
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length
            trimmed = seq[length:]
            #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
            #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant IDK
            yield {
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": direction,
                "length": seq_length-length
            }
    else:
        print(f"The length of your forward primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {forward_length}")




def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """