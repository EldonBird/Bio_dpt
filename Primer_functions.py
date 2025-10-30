import pandas as pd
import re 
from Bio.Seq import Seq
import logging
from typing import Dict
import primer3
from collections.abc import Callable
# from itertools import filter

# the order is 
# generate allele specific (and evaluate)
# generate matching (and evaluate)
# filter
# rank


# # Sequence -> Sequence
# def Reverse_Complement(sequence: str): #  -> str
    
#     return str(Seq(sequence).reverse_complement())
# This function was unnecessary


def introduce_mismatch(primer_sequence: str) -> str:
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
def calc_gc_content(sequence: str):
    gc_total = 0
    for nucleotide in sequence:
        if nucleotide == 'G' or nucleotide == 'C':
            gc_total += 1
    return gc_total/len(sequence)

def evaluate_primers(primer_dict: dict) -> Dict:
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
        gc_content = calc_gc_content(primer_dict['primer_sequence'])
        hairpin = primer3.bindings.calc_hairpin(primer_dict['primer_sequence'])
        homodimer = primer3.bindings.calc_homodimer(primer_dict['primer_sequence'])

        primer_dict["tm"] = tm_result
        primer_dict["gc_content"] = gc_content
        primer_dict["hairpin_dg"] = hairpin.dg 
        print(hairpin)
        primer_dict["homodimer_dg"] = homodimer.dg 
        print(homodimer)
        primer_dict["success"] = True
        return primer_dict
    except Exception as e:
        logging.error(f"Primer evaluation failed for sequence: {primer_dict['primer_sequence']}\nError: {str(e)}")
        primer_dict["tm"] = 0
        primer_dict["gc_content"] = 0
        primer_dict["hairpin_dg"] = 999
        primer_dict["homodimer_dg"] = 999
        primer_dict["success"] = False
        print(f"I have failed, please forgive my sins ({primer_dict['snpID']}, {primer_dict['allele']}, {primer_dict['direction']})")   
        return primer_dict

def rank_primers(primers: list[dict], target_tm = 62.5, target_gc = 50, optimism = 5) -> list[dict]:
    """
        Rank primers based on Tm proximity to 62.5Â°C and GC content.
        TODO: Refine ranking criteria.
        - Consider weighting Tm vs. GC scores.
        - Add user-configurable ranking metrics.
        """
    for primer in primers:
        primer["tm_score"] = abs(primer["tm"] - target_tm)
        primer["gc_score"] = abs(primer["gc_content"] - target_gc)
        primer["score"] = primer["tm_score"] + primer["hairpin_dg"] + primer["homodimer_dg"] + primer["gc_score"] 


    primers.sort(key=lambda x: x["score"])

    
    return primers[:optimism]
































# def metrics_for_list(primer_list: List[str], evaluate: Callable[[str], Dict[str, float]]) -> pd.DataFrame:
#     """
#     Build metrics table for a list of primers.
#     The columns are primer, tm, gc, hairprin, homodimer
#     In R, we were doing 
#             sapply(list_of_primers, calculate_tm)
#             sapply(list_of_primers, calculate_hairpin)
#             sapply(list_of_primers, calculate_homodimer)
#     then we indexed those results to apply the thresholds per list.
#     the calculate.... is used in the evaluate_primer function.           
#     """

#     if not primer_list:
#         return pd.DataFrame(columns=["primer", "tm", "gc", "hairpin", "homodimer"])
    
#     rows = []

#     for p in primer_list:
#         m = evaluate(p)
#         rows.append({
#             "primer": p,
#             "tm": float(m["tm"]),
#             "gc": float(m.get("gc_content", m.get("gc", float("nan")))),
#             "hairpin": float(m["hairpin"]),
#             "homodimer": float(m["homodimer"]),
#         })
#     return pd.DataFrame(rows)

# def _soft_keep(allele_list: list[dict], 
#                predicate: value_to_pass_somehow, 
#                keep_at_least: int,
#                closeness_key: pd.Series) -> list[dict]:
#     """
#     If suficient rows pass the predicate, then only keep those
#     Else keep the best keep_at_least rows by closest to the threshold

#     R code to check.
#     k = candidates[ sapply(candidates, calculate_homodimer)[2,] < Homodimer ]
#     if (length(k) > 5) keep k
#     else keep the ~5 closest to Homodimer
#     """
#     passed = []
#     for primer in allele_list:
#         if predicate == True:
#             passed.append(primer)

#     if len(passed) >= keep_at_least:
#         return passed
#     if not allele_list:
#         return allele_list
    
#     k = min(keep_at_least, len(allele_list))

#     return (allele_list.assign(_close=closeness_key.abs())
#               .sort_values("_close")
#               .head(k)
#               .drop(columns="_close"))

def filter_one_list_soft(allele_list: list[dict],
                         desired_tm: float = 60.0,
                         diff: float = 3.0,
                         homodimer_goal: float = 3.0,
                         hairpin_goal: float = 3.0,
                         keep_at_least: int = 5) -> list[dict]:
    
    """
    Soft filter a single candidate list such as the stage1_filter behavior
    Order is homodimer, hairpin, tm > lower, tm < upper
    
    Applies the 4 checks and then returns filtered list of primer string
    for that one candidate list

    The checks are homodimer, hairpin, tm, 

    Used in the R stage1_filter:
        homodimer < homodimer_goal
        hairpin   < hairpin_goal
        Tm        < desired_tm + diff      # "above upper" trim
        Tm        > desired_tm - diff      # "below lower" trim
    """
    low_fail = 0
    high_fail = 0
    dimer_fail = 0
    hairpin_fail = 0

    if not allele_list:
        return allele_list
    snp_allele = allele_list[0]['snpID'] + "_" + allele_list[0]['allele']


        # tm > min
    # allele_list_phair.sort(key=lambda x: x["tm"])
    allele_list_pltm = list(filter(lambda x : x["tm"] >= (desired_tm - diff), allele_list))
    if not allele_list_pltm:
        allele_list_pltm = allele_list
        print("No primers passed low bound Tm filter; proceeding with post hair!!!!!!!")
        for primer in allele_list:
            print({snp_allele}, {primer["tm"]})
        low_fail += 1
    # print(f"number of primers for {snp_allele} after lower bound Tm filter: {len(allele_list_pltm)}")

    # tm < max
    # allele_list_pltm.sort(reverse = True, key=lambda x: x["tm"])
    allele_list_phtm = list(filter(lambda x : x["tm"] <= (desired_tm + diff), allele_list_pltm))
    if not allele_list_phtm:
        allele_list_phtm = allele_list_pltm
        print("No primers passed upper bound Tm filter; proceeding with post lower bound tm!!!!!!!")
        for primer in allele_list_pltm:
            print({snp_allele}, {primer["tm"]})
        high_fail += 1
    # print(f"number of primers for {snp_allele} after upper bound Tm filter: {len(allele_list_phtm)}")

    
        # homodimer < max
    # allele_list.sort(key=lambda x: x["homodimer_dg"])
    allele_list_phomo = list(filter(lambda x : (homodimer_goal*-1) < x["homodimer_dg"] < homodimer_goal, allele_list_phtm))
    if not allele_list_phomo:
        allele_list_phomo = allele_list_phtm
        print("No primers passed homodimer filter; proceeding with allele list.")
        for primer in allele_list_phtm:
            print({snp_allele}, {primer["homodimer_dg"]})
        dimer_fail += 1
    # print(f"number of primers for {snp_allele} after homodimer filter: {len(allele_list_phomo)}")

    # hairpin < max
    # allele_list_phomo.sort(key=lambda x: x["hairpin_dg"])
    allele_list_phair = list(filter(lambda x : (hairpin_goal*-1) < x["hairpin_dg"] < hairpin_goal, allele_list_phomo))
    if not allele_list_phair:
        allele_list_phair = allele_list_phomo
        print("No primers passed hairpin filter; proceeding with post homo.")
        for primer in allele_list_phomo:
            print({snp_allele}, {primer["hairpin_dg"]})
        hairpin_fail += 1
    # print(f"number of primers for {snp_allele} after hairpin filter: {len(allele_list_phair)}")



    # print(low_fail, high_fail)
    return (low_fail, high_fail, dimer_fail, hairpin_fail)
    # return allele_list_phtm[:min(keep_at_least, len(allele_list_phtm))]











# def filter_primers(primers: list[dict],
#                    desired_tm: float = 64.0,
#                    diff: float = 3.0,
#                    hairpin_goal: float = 45.0,
#                    homodimer_goal: float = 45.0,
#                    keep_at_least: int = 5) -> pd.DataFrame:
#     """
#     Filter applied per snpID, allele, direction.
#     Keeps at least a few best available primers per group, 
#     unless a group had none to start with.
#     reurns a flat datafram with metrics so the downstream remain unchanged

#     We are filtering per each SNP group and keeps some best available 
#     Group candidates and collects sequences into lists.
#     (In R, each row already carried lists. 
#     We have candidates as rows, so we group to recreate per-group lists.)

#     """

#     """
#     filter primers took in a allele_list (but we're going to loop it on list of dicts for simplicity's sake.)
#     so I think I just need filter one list soft looped
#     and _soft keep perhaps (maybe we do things differently than archalie)
#     and I don't think we need metrics _for_list because the list's should already be made from generate_allele_specific_primers
    
#     """
#     if not primers:
#         return primers

#     # Group candidates per SNP/allele/direction
#     grouped = (primers.groupby(["snpID", "allele", "direction"])
#                       .agg({"primer_sequence": list})
#                       .reset_index())
#     #Isaiah note: shouldn't be necessary because they are already grouped into lists at that level

#     # Apply soft filter to each group's list
#     # Like stage1_filter() on each list in R, 
#     # same staged thresholds and bottlenecks fallbacks
#     def _apply(row):
#         return filter_one_list_soft(
#             row["primer_sequence"],
#             evaluate=evaluate_primer,
#             desired_tm=desired_tm,
#             diff=diff,
#             homodimer_goal=homodimer_goal,
#             hairpin_goal=hairpin_goal,
#             keep_at_least=keep_at_least
#         )
         

#     grouped["kept_sequences"] = grouped.apply(_apply, axis=1)
#     grouped = grouped[grouped["kept_sequences"].map(len) > 0]
    
#     # drop rows that ended empty, like farway empty
#     if grouped.empty:
#         return pd.DataFrame(columns=list(primers.columns) + ["tm","gc_content","hairpin","homodimer"])

#     # Explode back to rows and attach metrics (so rank_primers() still works)
#     # raank_primers will still receive a flat table with one primer per row 
#     # with tm, gc_content, hairpin, and homodimer
#     # the get_filter added substrings_count, maybe if we distinguish near/far
#     # implement it? If not keep it simple.
#     rows = []
#     for _, r in grouped.iterrows():
#         for seq in r["kept_sequences"]:
#             m = evaluate_primer(seq)
#             rows.append({
#                 "snpID": r["snpID"],
#                 "allele": r["allele"],
#                 "direction": r["direction"],
#                 "primer_sequence": seq,
#                 "length": len(seq),
#                 "tm": m["tm"],
#                 "gc_content": m["gc_content"],
#                 "hairpin": m["hairpin"],
#                 "homodimer": m["homodimer"],
#             })
#     return pd.DataFrame(rows)













def generate_allele_specific_primers(snps_list: list[dict], min_len: int = 18, max_len: int = 28) -> list[list[list[dict]]]:
    # make primers makes a dictionary for every length of one direction of an allele for a SNP.
    # those dictionaries are stored in a list, so a list for forward and a list for backward
    # those lists are stored in another list, one for each SNP. 
    # then those are stored in one big list, a list of all SNP given this session.
    # list_of_every_SNP[list_of_forward_and_reverse_directions[list_of_dictionaries[dictionary_of_particular_length]]]
    """
        Generate allele-specific primers (forward/reverse) ending at the SNP.
        TODO: Optimize for large SNP sets.
        - Use parallel processing (e.g., multiprocessing) for many SNPs.
        - Add validation for sequence length and SNP position.
        """
    all_primers = []
    min_len -= 2 #don't know why but 2 and 1 have to be removed from the inputs to get the desired lengths
    max_len -= 1

    for snp_dict in snps_list:
        #why pass this in seperately when it's already in the dict?
        snp_id = snp_dict["snpID"]
        allele = snp_dict["allele"]
        sequence = snp_dict["sequence"]
        snp_pos = snp_dict["position"]

        forward = sequence[snp_pos - max_len :snp_pos+1]#this gets the largest segment.   
        forward_mismatch = introduce_mismatch(forward)
    
        reverse = str(Seq(sequence[snp_pos:snp_pos+max_len+1]).reverse_complement()) #creates a Biopython sequence, gets the reverse complement, and converts is back to a string
        reverse_mismatch = introduce_mismatch(reverse)

        this_allele_primers = (make_primers(forward_mismatch, min_len, max_len, snp_id, allele))\
                            + (make_primers(reverse_mismatch, min_len, max_len, snp_id, allele, "reverse"))# this make one list of dictionaries for both flanking directions
        all_primers.append(this_allele_primers) #this adds this list to the larger list
    return all_primers # this will return a list of lists of dictionaries. Each allele is a list. 
#[[snp1 allele1 dictionaries],[snp1 allele2 dictionaries],[snp2 allele2 dictionaries],[snp2 allele2 dictionaries]] each snp and allele are on the same level.
   

def make_primers(seq, min_len, max_len, snp_id, allele, direction="forward") -> list[dict]: 
    seq_length = len(seq)
    primers = []
    if seq_length >= min_len:
        for length in range(max_len-min_len):#possible bug if the forward mismatch is smaller than the minimum length
            
            primary_primer= {
                "snpID": snp_id,
                "allele": allele,
                "primer_sequence": seq[length:], #this is the trimmed length
                "direction": direction,
                "length": seq_length-length 
            }
            primers.append(evaluate_primers(primary_primer))#Evaluate primers was built to handle one dictionary at a time. Adding it here saves time by avoiding an extra loop. 
            #maybe avoid calling the function all together by passing evaluate_primers directly into make_primers
    else:
        print(f"The length of your forward primer wasn't long enough. \nYou needed one at least {min_len} long and it ended up only being {seq_length}")
    return primers



def generate_matching_primers(snp_data: pd.DataFrame, allele_specific_primers: pd.DataFrame, min_dist: int = 100, max_dist: int = 500): # -> pd.DataFrame::
    """
        Generate matching primers for top 5 allele-specific primers.
        TODO: Optimize primer pairing.
        - Use primer3-py's designPrimers for more efficient pairing.
        - Add checks for primer pair compatibility (e.g., Tm difference < 5Â°C).
        """
