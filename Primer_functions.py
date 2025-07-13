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


def Evaluate_Primers(primer_seq: str) -> Dict:
    """
    Evaluates a single primer using the modern `primer3.design_primers` method
    with thermodynamic parameter path for accurate ΔG calculation.
    """
    # Standardize sequence to uppercase for consistency.
    seq = primer_seq.upper()
    # Validate input: return default values for empty or non-DNA sequences.
    if not seq or not set(seq) <= {'A', 'T', 'G', 'C'}:
        return {"tm": 0, "gc_content": 0, "hairpin_dg": 0, "homodimer_dg": 0}

    # Use the robust, all-in-one 'design_primers' engine for evaluation.
    result = primer3.design_primers(
        # Define the sequence to be evaluated.
        seq_args={
            'SEQUENCE_TEMPLATE': seq,
        },
        # Specify global parameters for the evaluation task.
        global_args={
            'PRIMER_TASK': 'check_primers', # Instructs primer3 to check a provided primer.
            'PRIMER_EXPLAIN_FLAG': 1, # Provides detailed output.
            # CRITICAL: Provide the path to thermodynamic parameters for accurate ΔG calculation.
            'PRIMER_THERMODYNAMIC_PARAMETERS_PATH': primer3.THERMODYNAMIC_PARAMETERS_PATH,
        },
        # Provide the primer sequence to be checked.
        primer_args={'PRIMER_SEQUENCE': seq}
    )
    
    # Extract results from the output dictionary.
    # The 'bindings' return ΔG in cal/mol, so we convert it to kcal/mol by dividing by 1000.
    return {
        'tm': result.get('PRIMER_LEFT_0_TM', 0.0),
        'gc_content': result.get('PRIMER_LEFT_0_GC_PERCENT', 0.0),
        'hairpin_dg': result.get('PRIMER_LEFT_0_HAIRPIN_DG', 0.0) / 1000.0,
        'homodimer_dg': result.get('PRIMER_LEFT_0_SELF_ANY_DG', 0.0) / 1000.0
    }

def Rank_Primers(primers: pd.DataFrame): # -> pd.DataFrame:
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    ...

def Filter_Primers(
    primers: pd.DataFrame,
    tm_range: Tuple[float, float] = (60.0, 65.0),
    gc_range: Tuple[float, float] = (40.0, 60.0),
    hairpin_dg_min: float = -9.0, # ΔG threshold for hairpins (less negative is better).
    homodimer_dg_min: float = -9.0, # ΔG threshold for homodimers.
    use_fallback: bool = True
) -> pd.DataFrame:
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
        metrics_df = primers['primer_sequence'].apply(evaluate_primer).apply(pd.Series)
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
        (primers_with_metrics['homodimer_dg'] > homodimer_dg_min)
    )
    # Apply the mask to get the subset of primers that passed.
    strict_results = primers_with_metrics[strict_filter]

    # If any primers passed the strict filter, add a 'filter_level' column and return them.
    if not strict_results.empty:
        return strict_results.assign(filter_level='strict')

    # --- Fallback Logic ---
    # If no primers passed strict and fallback is enabled, try again with relaxed criteria.
    if use_fallback:
        print("WARNING: No primers passed strict filtering. Applying relaxed criteria.")
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


    forward = sequence[center - max_len + 1 :center + 1]#this gets the largest segment.   
    forward_mismatch = Introduce_Mismatch(forward)
    forward_length = len(forward_mismatch)
   

    if forward_length >= min_len:
        for length in range(max_len-(min_len-1)):#possible bug if the forward missmatch is smaller than the minimum length
            trimmed = forward_mismatch[length:]
            this_allele_primers.append({
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
    reverse = str(Seq(sequence[center + 1:center+max_len-1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
    reverse_mismatch = Introduce_Mismatch(reverse)
    reverse_length = len(reverse_mismatch)
    

    if reverse_length >= min_len:
        for length in range(max_len-(min_len-1)):
            trimmed = forward_mismatch[length:]
            this_allele_primers.append({
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": trimmed,
                "direction": "reverse",
                "length": reverse_length-length
            })
            
    else:
        print(f"The length of your reverse primer wasn't long enough. \n You needed one at least {min_len} long and it ended up only being {reverse_length}")
    return pd.DataFrame(this_allele_primers)

def Generate_Matching_Primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """
