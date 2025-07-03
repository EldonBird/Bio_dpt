import pandas as pd
import primer3
import random


def add(A, B):
    return A + B

def Check_Multiplex_Compatibility(primer_pairs: pd.DataFrame, heterodimer_max: float = 50.0): # -> [{score: 32.2, combination: (P1, P2, P3, P4)}, {...}]:
    """
    Check primer sets for multiplex compatibility.
    TODO: Enhance for multiple SNPs.
    - Extend to check cross-SNP interactions (current checks within SNP).
    - Optimize for large primer sets using batch dimer calculations.
    """

    def _get_snp_ids(primer_pairs: pd.DataFrame) -> list:
        """
        Extract unique SNP IDs from primer pairs DataFrame.
        """
        # Assuming 'snpId' is a column in the primer_pairs DataFrame
        return primer_pairs['snpId'].unique().tolist()
    
    def _get_primer_pair_list(primer_pairs: pd.DataFrame) -> list:
        """
        Extract primer pairs from DataFrame.
        """
        # return [tuple(item) for item in primer_pairs[['Forward', 'Reverse']].values]
        # TODO: add index from 0 to len(primer_pairs)
        return []
    
    def _calculate_hetero_compatibility_score(heterodimer_results) -> float:
        """
        Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
        """
        # Example scoring logic (to be replaced with actual logic)
        # score = 0.0
        # for result in hairpin_results + homodimer_results + heterodimer_results:
        #     if result['tm'] < heterodimer_max:
        #         score += 1.0
        # return score / len(hairpin_results + homodimer_results + heterodimer_results)
        return random.uniform(0, 100)  # Placeholder for actual score calculation
    
    def _calculate_homo_compatibility_score(hairpin_results, homodimer_results, heterodimer_results) -> float:
        """
        Calculate compatibility score based on hairpin, homodimer, and heterodimer results.
        """
        # Example scoring logic (to be replaced with actual logic)
        # score = 0.0
        # for result in hairpin_results + homodimer_results + heterodimer_results:
        #     if result['tm'] < heterodimer_max:
        #         score += 1.0
        # return score / len(hairpin_results + homodimer_results + heterodimer_results)
        return random.uniform(0, 100)  # Placeholder for actual score calculation

    ### 1. Prepare the data
    pp_list = _get_primer_pair_list(primer_pairs) # [(primer_forward, primer_reverse), ...]

    ### 2. Define compatibility score calculation for all primer pairs
    homo_table = [0] * len(pp_list) # score of homodimer&hairpin. as score goes higher, it makes worse effect on pcr
    hetero_table = [[0] * len(pp_list)] * len(pp_list) # score of heterodimer. as score goes higher, it makes worse effect on pcr

    for i in range(len(pp_list)):
        primer1_f, primer1_r = pp_list[i]

        # calculate hairpin&homodimer score
        hairpin_results = []
        hairpin_results.append(primer3.calcHairpin(primer1_f.sequence))
        hairpin_results.append(primer3.calcHairpin(primer1_r.sequence))
        hairpin_results.append(primer3.calcHairpin(primer2_f.sequence))
        hairpin_results.append(primer3.calcHairpin(primer2_r.sequence))

        homodimer_results = []
        homodimer_results.append(primer3.calcHomodimer(primer1_f.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer1_r.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer2_f.sequence))
        homodimer_results.append(primer3.calcHomodimer(primer2_r.sequence))

        homo_table[i] = _calculate_homo_compatibility_score(hairpin_results, homodimer_results)

        for j in range(i + 1, len(pp_list)):
            if(i == j):
                continue
            primer2_f, primer2_r = pp_list[j]
            # Calculate heterodimer score
            heterodimer_results = []
            heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer1_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_f.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_f.sequence, primer2_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_f.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer1_r.sequence, primer2_r.sequence))
            heterodimer_results.append(primer3.calcHeterodimer(primer2_f.sequence, primer2_r.sequence))

            result = _calculate_hetero_compatibility_score(heterodimer_results)
            hetero_table[i][j] = result
            hetero_table[j][i] = result

    ### 3. calculate compatibility scores for all **combinations** of primer pairs
    # prepare the all combinations of primer pairs
    from itertools import product
    snp_ids = _get_snp_ids(primer_pairs)  # Extract unique SNP IDs from primer pairs
    grouped = []
    for snp in snp_ids:
        grouped.append([pp for pp in pp_list if pp['snpId'] == snp])

    combinations = list(product(*grouped)) # [(P1, P2, P3, P4), (P5, P2, P3, P4)] if len(snp_ids) == 4 

    # calculate compatibility scores for all combinations
    scores = []
    for comb in combinations: # comb is (P1, P2, P3, P4)
        score = 0
        for i in range(len(comb)):
            pp1 = comb[i]
            score += homo_table[pp1.index]
            for j in range(i + 1, len(comb)):
                pp2 = comb[j]
                score += hetero_table[pp1.index][pp2.index]
        scores.append(score)
    
    # sort combinations by their score
    zipped = list(zip(scores, combinations))
    zipped_sorted = sorted(zipped, key=lambda x: x[0])
    sorted_scores, sorted_combinations = zip(*zipped_sorted)
    
    ### Packing result
    result = []
    for i in range(len(sorted_scores)):
        dict = {}
        dict['score'] = sorted_scores[i]
        dict['combination'] = sorted_combinations[i]
        result.append(dict)
    return result

    