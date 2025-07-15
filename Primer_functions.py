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


def Evaluate_Primers(primer_seq: str): # -> Dict:
    """
        Evaluate primer quality using primer3-py.
        TODO: Enhance error handling.
        - Handle primer3-py failures gracefully.
        - Add logging for failed evaluations.
        """
    # print(f'this is the primer given {primer_seq}')
    # # Setup logging
    # logging.basicConfig(
    #     # print('in set up logic')
    #     filename="primer_evaluation.log",
    #     level=logging.INFO,
    #     # print('logging info')
    #     format="%(asctime)s - %(levelname)s - %(message)s"
    # )

    # try:
    #     result = primer3.bindings.design_primers(
    #         {
    #             "SEQUENCE_TEMPLATE": primer_seq,
    #             "SEQUENCE_PRIMER": primer_seq
    #         },
    #         {
    #             "PRIMER_OPT_SIZE": 20,
    #             "PRIMER_MIN_SIZE": 18,
    #             "PRIMER_MAX_SIZE": 28,
    #             "PRIMER_OPT_TM": 62.5,
    #             "PRIMER_MIN_TM": 60.0,
    #             "PRIMER_MAX_TM": 65.0,
    #             "PRIMER_MIN_GC": 40.0,
    #             "PRIMER_MAX_GC": 60.0,
    #             "PRIMER_MAX_HAIRPIN_TH": 45.0,
    #             "PRIMER_MAX_SELF_ANY_TH": 45.0
    #         }
    #     )
    #     return {
    #         "tm": result.get("PRIMER_LEFT_0_TM", 0),
    #         "gc_content": result.get("PRIMER_LEFT_0_GC_PERCENT", 0),
    #         "hairpin": result.get("PRIMER_LEFT_0_HAIRPIN_TH", 999),
    #         "homodimer": result.get("PRIMER_LEFT_0_SELF_ANY_TH", 999),
    #         "success": True
    #     }
    # except Exception as e:
    #     logging.error(f"Primer evaluation failed for sequence: {primer_seq}\nError: {str(e)}")
    #     return {
    #         "tm": 0,
    #         "gc_content": 0,
    #         "hairpin": 999,
    #         "homodimer": 999,
    #         "success": False
    #     }
        # print(f'this is the primer given {primer_seq}')
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
        result = primer3.bindings.calc_tm(primer_seq)
        # gc_content = primer3.bindings.calc_gc(primer_seq)
        hairpin = primer3.bindings.calc_hairpin(primer_seq)
        homodimer = primer3.bindings.calc_homodimer(primer_seq)
        
        return {
            "tm": result,
            # "gc_content": gc_content,
            "hairpin": hairpin.tm if hasattr(hairpin, 'tm') else 0,
            "homodimer": homodimer.tm if hasattr(homodimer, 'tm') else 0,
            "success": True
        }
    except Exception as e:
        logging.error(f"Primer evaluation failed for sequence: {primer_seq}\nError: {str(e)}")
        print(primer_seq)
        return {
            "tm": 0,
            "gc_content": 0,
            "hairpin": 999,
            "homodimer": 999,
            "success": False
        }

def rank_primers(primers: pd.DataFrame, target_tm = 62.5, target_gc = 50) -> pd.DataFrame:
    """
    Rank primers based on Tm proximity to 62.5Â°C and GC content.
    TODO: Refine ranking criteria.
    - Consider weighting Tm vs. GC scores.
    - Add user-configurable ranking metrics.
    """
    primers["tm_score"] = abs(primers["tm"] - target_tm)
    primers["gc_score"] = abs(primers["gc_content"] - target_gc)
    primers["score"] = primers["tm_score"] + primers["gc_score"] + primers["hairpin"] + primers["homodimer"]
    return primers.sort_values("score").groupby(["snpID", "allele", "direction"]).head(5)
    Rank primers based on Tm proximity to 62.5Â°C and GC content.
    TODO: Refine ranking criteria.
    - Consider weighting Tm vs. GC scores.
    - Add user-configurable ranking metrics.
    """
    primers["tm_score"] = abs(primers["tm"] - target_tm)
    primers["gc_score"] = abs(primers["gc_content"] - target_gc)
    primers["score"] = primers["tm_score"] + primers["gc_score"] + primers["hairpin"] + primers["homodimer"]
    return primers.sort_values("score").groupby(["snpID", "allele", "direction"]).head(5)
    ...

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
    #made these two separate functions incase we want to parallelize them
    for snp_dict in snps_list:
        all_primers.append(Find_Primers(snp_dict, min_len, max_len))#call the function again and again. My say is we batch by snpID or something and multiprocess a batch
        #starting a whole thread just for a couple for loops doesn't quite seem justified.
    return all_primers
    # make an empty list of primers
    #iterate through each function calling the find primers function. 
    #this function is a paralizing shell that will house the true find primers function.
    return all_primers

def Find_Primers(snp, min_len, max_len):
    this_allele_primers = []#a list of dictionaries
    snp_id = snp["snpID"]
    allele = snp["allele"]
    sequence = snp["sequence"]
    snp_pos = snp["position"]

    # print(f'the snp: {sequence[(snp_pos-3):(snp_pos+5)]}')
    # print(f'the snp Gone?: {sequence[:snp_pos]}')
    print(f'the snp alone: {sequence[snp_pos]}')
    forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
    forward_mismatch = Introduce_Mismatch(forward)
    forward_length = len(forward_mismatch)
    # print(f'the forward: {forward}')

    if forward_length >= min_len:
        for length in range(max_len-min_len+1):#possible bug if the forward missmatch is smaller than the minimum length
            # print(forward_mismatch[length:])
            trimmed = forward_mismatch[length:]
            # print(f'here\'s the sequence: {trimmed}')
            this_allele_primers.append({#take this part out of the loop, so we can have one dictioanry that says the SNP ID and ALLELE and Dirrection, and then a list in that 
                #dictioanry of sequenec and lengths. Storing the name over and over seems redundednt IDK
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": "forward",
                "length": forward_length-length
            })
            
    else:
        print(f"The length of your forward primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {forward_length}")


    # Reverse primer: downstream sequence, reverse complemented
    # print(f"here's the sequence again: {sequence}")
    # print(f'the thing before modify: {sequence[snp_pos:snp_pos+max_len+1]}')
    reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
    reverse_mismatch = Introduce_Mismatch(reverse)

    reverse_length = len(reverse_mismatch)
    
    # print(f'the reverse: {reverse}')
    if reverse_length >= min_len:
        for length in range(max_len-min_len+1):
            trimmed = reverse_mismatch[length:]
            this_allele_primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": "reverse",
                "length": reverse_length-length
            })

            
    else:
        print(f"The length of your reverse primer wasn't long enough. \n You needed one at least {min_len} long and it ended up only being {reverse_length}")
    return this_allele_primers

def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """