import pandas as pd
import re 
from Bio.Seq import Seq
import logging
from typing import Dict
import primer3

# the order is 
# generate allele specific (and evaluate)
# generate matching (and evaluate)
# filter
# rank


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
        print("Warning: Invalid primer input.") #Correction needed: Use logging instead of print for consistency.
        return primer_sequence

    primer_sequence = primer_sequence.upper().strip()

    # Must only contain A, C, G, T
    if not re.match("^[ACGT]+$", primer_sequence):
        print(f"Warning: Invalid characters in primer: {primer_sequence}") #Correction needed: Use logging.
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
def Calc_GC_Content(sequence: str): #Could use BioPython's GC function for efficiency
    gc_total = 0
    for nucleotide in sequence:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_total += 1
    return gc_total/len(sequence) #Correction needed: Handle empty sequence to avoid division by zero.

def Evaluate_Primers(primer_dict: dict) -> Dict:
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
        gc_content = Calc_GC_Content(primer_dict['primer_sequence'])
        hairpin = primer3.bindings.calc_hairpin(primer_dict['primer_sequence'])
        homodimer = primer3.bindings.calc_homodimer(primer_dict['primer_sequence'])

        primer_dict["tm"] = tm_result
        primer_dict["gc_content"] = gc_content
        primer_dict["hairpin_dg"] = hairpin.tm if hasattr(hairpin, 'tm') else 0 
        # ^Correction needed: Use dg instead of tm for accuracy; primer3 hairpin returns structure with dg.
        primer_dict["homodimer_dg"] = homodimer.tm if hasattr(homodimer, 'tm') else 0
        # ^Correction needed: Same as above, use dg.
        primer_dict["success"] = True
        return primer_dict
    except Exception as e:
        logging.error(f"Primer evaluation failed for sequence: {primer_dict['primer_sequence']}\nError: {str(e)}")
        primer_dict["tm"] = 0
        primer_dict["gc_content"] = 0
        primer_dict["hairpin_dg"] = 999
        primer_dict["homodimer_dg"] = 999
        primer_dict["success"] = False
        return primer_dict

def rank_primers(primers: list[dict], target_tm = 62.5, target_gc = 50, optimism = 5) -> list[dict]:
    """
        Rank primers based on Tm proximity to 62.5°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """ # ^Correction Needed: add heterodimer to score.
    for primer in primers:
        primer["tm_score"] = abs(primer["tm"] - target_tm)
        primer["gc_score"] = abs(primer["gc_content"] - target_gc)
        primer["score"] = primer["tm_score"] + primer["hairpin_dg"] + primer["homodimer_dg"] + primer["gc_score"] 
        # ^Correction needed: Add heterodimer_dg; use negative DG for penalties (less negative is better).


    primers.sort(key=lambda x: x["score"])

    
    return primers[:optimism]

def Filter_Primers(
        #Correction Needed: Correction needed: Switch to list[dict] input/output as per comments; implement per-group filtering.
    # This is not filtering right. It's filtering the whole dataframe, and then if nothing passes, it's relaxing the requirements and trying again.
    # We need it to filter per list of dictionaries (per allele flanking direction). 
    # It should filter the primers close to the SNP and if none of them work then just take the best 5. 
    # Then the ones far away and if none of them work, scrap the whole allele flanking direction.
    # We also need it to use list of dictionaries rather than pandas dataframes. 
    # https://claude.ai/share/e73c5d8d-0209-485f-aff6-9db0f6066dd4 My claude conversation on the topic
    primers: pd.DataFrame, #Correction needed: Change to list[list[list[dict]]].
    tm_range: tuple[float, float] = (60.0, 65.0),
    gc_range: tuple[float, float] = (40.0, 60.0),
    hairpin_dg_min: float = -9.0, # ΔG threshold for hairpins (less negative is better).
    homodimer_dg_min: float = -9.0, # ΔG threshold for homodimers.
    use_fallback: bool = True
) -> pd.DataFrame: #Correction needed: Change to list.
    """
    Evaluates (if necessary) and filters primers based on quality metrics.
    - Tm must be within tm_range.
    - GC% must be within gc_range.
    - Hairpin and homodimer ΔG must be weaker than the threshold (less negative, i.e., > min).
    """
    # Immediately return if the input DataFrame is empty.
    if primers.empty:
        return pd.DataFrame()
    
    # --- Evaluation Step ---
    # Define the set of required metric columns.
    required_cols = {'tm', 'gc_content', 'hairpin_dg', 'homodimer_dg'}
    if not required_cols.issubset(primers.columns):
        print("Metrics not found, running evaluation...")
        # Calculate metrics for each primer sequence.
        metrics_df = primers['primer_sequence'].apply(Evaluate_Primers).apply(pd.Series) #Correction needed: Add heterodimer calc.
        # Join the new metrics back to the original primer data.
        primers_with_metrics = primers.join(metrics_df)
    else:
        # If metrics are already present, just use the input DataFrame.
        primers_with_metrics = primers

    # --- Strict Filtering Logic ---
    # Create a boolean mask where each condition must be True for a primer to pass.
    strict_filter = (
        primers_with_metrics['tm'].between(*tm_range) &
        primers_with_metrics['gc_content'].between(*gc_range) &
        (primers_with_metrics['hairpin_dg'] > hairpin_dg_min) &
        (primers_with_metrics['homodimer_dg'] > homodimer_dg_min) #Correction needed: Add heterodimer check.
    )
    # Apply the mask to get the subset of primers that passed.
    strict_results = primers_with_metrics[strict_filter]

    # If any primers passed the strict filter, add a 'filter_level' column and return them.
    if not strict_results.empty:
        return strict_results.assign(filter_level='strict')

    # --- Fallback Logic ---
    # If no primers passed strict and fallback is enabled, try again with relaxed criteria.
    if use_fallback:
        print("WARNING: No primers passed strict filtering. Applying relaxed criteria.") #Correction needed: Use logging; make per-group.
        # Create a new boolean mask with wider, more tolerant thresholds.
        relaxed_filter = (
            primers_with_metrics['tm'].between(tm_range[0] - 2.0, tm_range[1] + 2.0) &
            primers_with_metrics['gc_content'].between(gc_range[0] - 5.0, gc_range[1] + 5.0) &
            (primers_with_metrics['hairpin_dg'] > hairpin_dg_min - 2.0) &
            (primers_with_metrics['homodimer_dg'] > homodimer_dg_min - 2.0)
        )
        # Apply the relaxed filter.
        relaxed_results = primers_with_metrics[relaxed_filter]

        # If any primers passed the relaxed filter, return them.
        if not relaxed_results.empty:
            return relaxed_results.assign(filter_level='relaxed')

    # If no primers pass even the relaxed criteria, print a warning and return an empty DataFrame.
    print("WARNING: No primers passed filtering, even with relaxed criteria.")
    return pd.DataFrame()


def Generate_Allele_Specific_Primers(snps_list: list[dict], min_len: int = 18, max_len: int = 28) -> list[list[list[dict]]]:
    # make primers makes a dictionary for every length of one direction of an SNP.
    # those dictionaries are stored in a list, so a list for forward and a list for backward
    # those lists are stored in another list, one for each SNP. 
    # then those are stored in one big list, a list of all SNP given this session.
    # list_of_every_SNP[list_of_forward_and_reverse_directions[list_of_dictionaries[dictionary_of_particular_length]]]
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        - Add restriction sites
        """
    all_primers = []
    min_len -= 2 #don't know why but 2 and 1 have to be removed from the inputs to get the desired lengths
    # ^ Correction needed: Remove arbitrary adjustment; explain or fix bug.
    max_len -= 1 #Correction needed: Same as above.

    for snp_dict in snps_list:
        this_allele_primers = []#a list of dictionaries
        snp_id = snp_dict["snpID"]
        allele = snp_dict["allele"]
        sequence = snp_dict["sequence"]
        snp_pos = snp_dict["position"]

        forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
        #^ Correction needed: Add length validation.
        forward_mismatch = Introduce_Mismatch(forward)
    
        reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
        reverse_mismatch = Introduce_Mismatch(reverse)

        this_allele_primers.append(Make_Primers(forward_mismatch, min_len, max_len, snp_id, allele))
        this_allele_primers.append(Make_Primers(reverse_mismatch, min_len, max_len, snp_id, allele, "reverse"))
 
        all_primers.append(this_allele_primers)
    return all_primers
   

def Make_Primers(seq, min_len, max_len, snp_id, allele, direction="forward") -> list[dict]: 
    seq_length = len(seq)
    primers = []
    if seq_length >= min_len:
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length
            #Incorrect: Bug if seq < min_len after mismatch; incorrect range (excludes max_len). Correction: Use range(min_len, max_len+1); validate seq length post-mismatch.
            trimmed = seq[length:] #Incorrect: No reverse-specific trimming logic. Correction: Trim differently for reverse like seq[:length]).
            #take this part out of the loop, so we can have one dictionary that says the SNP ID and ALLELE and Direction, 
            #and then a list in that dictionary of sequence and lengths. Storing the name over and over seems redundant IDK
            #Correction needed: For reverse, trim differently if needed.
            
            primary_primer= {
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": direction,
                "length": seq_length-length 
            } #Incorrect: Redundant fields; missing annealing positions, restriction sites. Correction: Add annealing_start/end, restriction_site keys.
            primers.append(Evaluate_Primers(primary_primer))#Evaluate primers was built to handle one dictionary at a time. Adding it here saves time by avoiding an extra loop. 
            #maybe avoid calling the function all together by passing Evaluate_Primers directly into make_primers
    else:
        print(f"The length of your forward primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {seq_length}") #Correction needed: Use logging; generalize for reverse.
    return primers



def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame:: 
    # Incorrect: Completely stubbed out; DataFrame-based inputs mismatch list outputs; no pair compatibility checks. Correction: Implement logic using primer3-py designPrimers, list inputs, Tm difference <5°C, heterodimer ΔG >-7 kcal/mol.
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5°C).
        """
    # Correction: implement logic to generate 5 pairs per top allele-specific primer, checking Tm and heterodimer.
# Correction: Add logic to extract sequences 100–500 bp away (upstream for reverse, downstream for forward), use primer3-py to design primers, and return list of pairs with metadata (Tm, GC, ΔG, annealing positions).