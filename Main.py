from Multiplex import *
from Output import *
from Primer_functions import *
import primer3
import re # run 'pip install regex' if not already installed
import time # to handle rate limiting
import requests
from typing import List # for type hinting

# Ensembl REST API base URL
ENSEMBL_REST = "https://rest.ensembl.org"

def Fetch_SNP_Data(rsids: List[str], flank_length: int = 800) -> list[dict]:
    """
    Retrieves SNP data from the Ensembl API, including flanking sequences and alleles.

    Args:
        rsids: List of SNP identifiers (e.g., ["rs1799971"]).
        flank_length: Number of base pairs to include on either side of the SNP.

    Returns:
        DataFrame containing SNP ID, allele, modified sequence, and SNP position.
    """
    snp_data = []
    headers = {"Content-Type": "application/json"}  # Required by Ensembl for JSON responses

    for rsid in rsids:
        try:
            # Step 1: Get SNP mapping (chromosome + position info)
            var_resp = requests.get(f"{ENSEMBL_REST}/variation/homo_sapiens/{rsid}?", headers=headers)
            print(f"getting info on {rsid}")
            var_resp.raise_for_status()
            var_data = var_resp.json()

            mappings = var_data.get("mappings", [])
            if not mappings:
                print(f"Warning: No mappings found for {rsid}.")
                continue

            # We'll just take the first mapping (usually sufficient for common SNPs)
            mapping = mappings[0]
            chrom = mapping["seq_region_name"]  # e.g., "11"
            pos = int(mapping["start"])   

            # Extract allele string like "A/G" or "C/T"
            allele_str = mapping.get("allele_string", "")
            ancestral = mapping.get("ancestral_allele")

            #Isaiah change: I dropped any ancestral Alleles that we knew were normal
            #In the future we should pull the list, split it, and have the user pick which ones they want.
            alleles = allele_str.split("/") if allele_str else []
            if ancestral in ['A', 'C', 'G', 'T']:
                alleles.remove(ancestral)


            # Ensure there are at least two alleles to work with
            if len(alleles) < 1:
                print(f"Warning: Less than 2 alleles for {rsid}. Skipping.")
                continue

            # Step 2: Fetch the flanking DNA sequence around the SNP
            seq_start = max(1, pos - flank_length)  # 1-based for Ensemble 
            seq_end = pos + flank_length #might run off the end of the chromosome if very unlucky
            seq_url = f"{ENSEMBL_REST}/sequence/region/human/{chrom}:{seq_start}..{seq_end}:1?"
            # print (f'chrom: {chrom}\nseq_start: {seq_start}\nseq_end: {seq_end}')
            seq_resp = requests.get(seq_url, headers={"Content-Type": "text/plain"})
            
            seq_resp.raise_for_status()
          
            template_seq = seq_resp.text.strip()

            # Position of the SNP relative to the start of the fetched sequence
            rel_pos = flank_length if seq_start > 1 else pos #should just be flanking length

            # Step 3: Replace the SNP base with each possible allele to simulate variation
            for allele in alleles:
                # Validate that allele contains valid DNA characters only
                if not re.fullmatch("[ACGTNacgtn]+", allele):
                    print(f"Skipping non-standard allele '{allele}' for {rsid}")
                    continue

                # Insert allele at the SNP site
                modified_seq = template_seq[:rel_pos] + allele.upper() + template_seq[rel_pos + 1:]

                # Append to results
                snp_data.append({
                    "snpID": rsid,
                    "allele": allele.upper(),
                    "sequence": modified_seq,
                    "position": rel_pos
                })

            # Sleep to respect Ensembl's rate limit (max 15 req/sec)
            time.sleep(0.34)

        except Exception as e:
            print(f"Error processing {rsid}: {e}")

    # If nothing was successfully retrieved, return an empty DataFrame
    if not snp_data:
        print("No valid SNP data could be retrieved.")
        return []
  

    return snp_data




def Main():
    """
        Main function to generate, filter, rank, pair, and export primers.
        TODO: Add comprehensive error handling and logging.
        - Log progress and errors to a file.
        - Add input validation for rsIDs and output format.
        TODO: Test with real SNP data.
        - Validate output with biological experts.
        - Benchmark performance for large SNP sets.
        """
    
    # real fetch snp
    # snp_df = Fetch_SNP_Data(["rs1799971", "rs12184297", "rs116801199", "rs12565286", "rs2977670", "rs28454925"], 30)# just here for testing.  , "rs599839"


    snp_df = [{'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, 
                {'snpID': 'rs12184297', 'allele': 'T', 'sequence': 'CTTTAAACCTCAACACATTATCAAGCATAATACTGTATATAATAAGTACTCAATACTGAAT', 'position': 30}, 
                {'snpID': 'rs116801199', 'allele': 'G', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATGAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
                {'snpID': 'rs116801199', 'allele': 'T', 'sequence': 'TAAAAAATGAATCTAATAATGAGGAAACATTAGAAAAAACCAAACTGAGGGATATTCTACA', 'position': 30}, 
                {'snpID': 'rs12565286', 'allele': 'G', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGGCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
                {'snpID': 'rs12565286', 'allele': 'C', 'sequence': 'GGAAGCATCCTTCACTATCTTCTACCAAGGCCTTCCTCCTTTGGTGCTTCAAAATTTTTTA', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'G', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGGGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'A', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGAGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'C', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGCGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs2977670', 'allele': 'T', 'sequence': 'AACCTTGGAGGACCTATTGCTTAAGGTGTGTGCCAAAGAAAGTAAGTTAGGGCAAGAGACT', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'C', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTCGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'G', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTGGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}, 
                {'snpID': 'rs28454925', 'allele': 'T', 'sequence': 'GGATTCGAATGGAAAGACATGGAATGGACTTGATTGGAATGGGTTGGGATGGAATGATCTA', 'position': 30}]
    # print(snp_df)

    primers = generate_allele_specific_primers(snp_df, 24, 30)
    # for prime_list in primers:
    #     for primer in prime_list:
    #         print(primer)
    #     print()
    # print(f"number of snps {len(primers)}")
    # print(f"number of directions {len(primers[0])}")
    # print(f"number of primers in a direction {len(primers[0][0])}")
    # print(len(snp_df))
    low = 0
    high = 0
    dimer = 0
    hairpin = 0
    for allele in primers:
        
        low_fail, high_fail, dimer_fail, hairpin_fail = filter_one_list_soft(allele, diff = 5.0)
        low += low_fail
        high += high_fail
        dimer += dimer_fail
        hairpin += hairpin_fail
    print(low)
    print(high)
    print(dimer)
    print(hairpin)

if(__name__ == "__main__"):
    Main()

# [[{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 26, 'tm': 62.64091186246884, 'gc_content': 0.46153846153846156, 'hairpin_dg': 40.246484286148586, 'homodimer_dg': 11.618441795495414, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 25, 'tm': 60.741963776567445, 'gc_content': 0.44, 'hairpin_dg': 39.085883205265475, 'homodimer_dg': -11.318504829494259, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 24, 'tm': 59.65329123118744, 'gc_content': 0.4166666666666667, 'hairpin_dg': 36.03872822148651, 'homodimer_dg': -32.43413848254389, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 26, 'tm': 69.86363516091336, 'gc_content': 0.5769230769230769, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 27.76279186528285, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 25, 'tm': 69.765715142391, 'gc_content': 0.6, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 25.839909279650783, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 24, 'tm': 68.83426028792064, 'gc_content': 0.625, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 24.30006522953198, 'success': True}], 
# [{'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 26, 'tm': 63.41694547346191, 'gc_content': 0.5, 'hairpin_dg': 39.73007908831477, 'homodimer_dg': -11.661466527887058, 'success': True}, 
#  {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 25, 'tm': 61.58199577952456, 'gc_content': 0.48, 'hairpin_dg': 33.026282568142335, 'homodimer_dg': -22.267119198164437, 'success': True}, 
#  {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 24, 'tm': 60.560716263151505, 'gc_content': 0.4583333333333333, 'hairpin_dg': 0.0, 'homodimer_dg': -30.140570030495894, 'success': True}, 
#  {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 26, 'tm': 70.08218467351901, 'gc_content': 0.6153846153846154, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 27.76279186528285, 'success': True}, 
#  {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 25, 'tm': 69.99393494968808, 'gc_content': 0.64, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 25.839909279650783, 'success': True}, 
# {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 24, 'tm': 69.08012732025867, 'gc_content': 0.6666666666666666, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 24.30006522953198, 'success': True}]]






# [[{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 26, 'tm': 62.64091186246884, 'gc_content': 0.46153846153846156, 'hairpin_dg': 40.246484286148586, 'homodimer_dg': 11.618441795495414, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 25, 'tm': 60.741963776567445, 'gc_content': 0.44, 'hairpin_dg': 39.085883205265475, 'homodimer_dg': -11.318504829494259, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 24, 'tm': 59.65329123118744, 'gc_content': 0.4166666666666667, 'hairpin_dg': 36.03872822148651, 'homodimer_dg': -32.43413848254389, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 26, 'tm': 69.86363516091336, 'gc_content': 0.5769230769230769, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 27.76279186528285, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 25, 'tm': 69.765715142391, 'gc_content': 0.6, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 25.839909279650783, 'success': True}, 
#   {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 24, 'tm': 68.83426028792064, 'gc_content': 0.625, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 24.30006522953198, 'success': True}], 
#   [{'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 26, 'tm': 63.41694547346191, 'gc_content': 0.5, 'hairpin_dg': 39.73007908831477, 'homodimer_dg': -11.661466527887058, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 25, 'tm': 61.58199577952456, 'gc_content': 0.48, 'hairpin_dg': 33.026282568142335, 'homodimer_dg': -22.267119198164437, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACG', 'direction': 'forward', 'length': 24, 'tm': 60.560716263151505, 'gc_content': 0.4583333333333333, 'hairpin_dg': 0.0, 'homodimer_dg': -30.140570030495894, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 26, 'tm': 70.08218467351901, 'gc_content': 0.6153846153846154, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 27.76279186528285, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 25, 'tm': 69.99393494968808, 'gc_content': 0.64, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 25.839909279650783, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'G', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATC', 'direction': 'reverse', 'length': 24, 'tm': 69.08012732025867, 'gc_content': 0.6666666666666666, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 24.30006522953198, 'success': True}], 
#    [{'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'AACCTCAACACATTATCAAGCATGAC', 'direction': 'forward', 'length': 26, 'tm': 60.3496111205028, 'gc_content': 0.38461538461538464, 'hairpin_dg': 38.69037508338931, 'homodimer_dg': -9.158301362534303, 'success': True}, 
#     {'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'ACCTCAACACATTATCAAGCATGAC', 'direction': 'forward', 'length': 25, 'tm': 59.87554310584028, 'gc_content': 0.4, 'hairpin_dg': 38.69037508338931, 'homodimer_dg': -9.158301362534303, 'success': True}, 
#     {'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'CCTCAACACATTATCAAGCATGAC', 'direction': 'forward', 'length': 24, 'tm': 58.42404345889554, 'gc_content': 0.4166666666666667, 'hairpin_dg': 38.69037508338931, 'homodimer_dg': -9.158301362534303, 'success': True}, 
#     {'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'GTATTGAGTACTTATTATATACAATG', 'direction': 'reverse', 'length': 26, 'tm': 49.94236197513237, 'gc_content': 0.23076923076923078, 'hairpin_dg': 0.0, 'homodimer_dg': -7.8920666413353615, 'success': True}, 
#     {'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'TATTGAGTACTTATTATATACAATG', 'direction': 'reverse', 'length': 25, 'tm': 48.32824792287022, 'gc_content': 0.2, 'hairpin_dg': 0.0, 'homodimer_dg': -7.8920666413353615, 'success': True}, 
#     {'snpID': 'rs12184297', 'allele': 'C', 'primer_sequence': 'ATTGAGTACTTATTATATACAATG', 'direction': 'reverse', 'length': 24, 'tm': 48.07499974324219, 'gc_content': 0.20833333333333334, 'hairpin_dg': 0.0, 'homodimer_dg': -7.8920666413353615, 'success': True}], 


{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CCTGGGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 30, 'tm': 67.57144768615979, 'gc_content': 0.5, 'hairpin_dg': -818.1391446079033, 'homodimer_dg': -6077.7359759805695, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CTGGGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 29, 'tm': 66.07258913299933, 'gc_content': 0.4827586206896552, 'hairpin_dg': -818.1391446079033, 'homodimer_dg': -6077.7359759805695, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TGGGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 28, 'tm': 65.58368362782704, 'gc_content': 0.4642857142857143, 'hairpin_dg': -818.1391446079033, 'homodimer_dg': -6077.7359759805695, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 27, 'tm': 64.41334458717705, 'gc_content': 0.48148148148148145, 'hairpin_dg': -703.3579156894302, 'homodimer_dg': -5256.433518143604, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 26, 'tm': 62.64091186246884, 'gc_content': 0.46153846153846156, 'hairpin_dg': -431.97164461101056, 'homodimer_dg': -2748.2905181436217, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 25, 'tm': 60.741963776567445, 'gc_content': 0.44, 'hairpin_dg': -291.93664461100707, 'homodimer_dg': -3446.5361446110037, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 24, 'tm': 59.65329123118744, 'gc_content': 0.4166666666666667, 'hairpin_dg': 81.7670421568364, 'homodimer_dg': -2504.5436867678472, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GCGGTTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 30, 'tm': 74.83615011034732, 'gc_content': 0.6333333333333333, 'hairpin_dg': -4695.662915689434, 'homodimer_dg': -14007.483518143592, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGTTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 29, 'tm': 73.25391354307897, 'gc_content': 0.6206896551724138, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -8728.631060300439, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGTTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 28, 'tm': 71.80362762506468, 'gc_content': 0.6071428571428571, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -7452.2014338237495, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GTTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 27, 'tm': 70.3716785611486, 'gc_content': 0.5925925925925926, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -7452.2014338237495, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 26, 'tm': 69.86363516091336, 'gc_content': 0.5769230769230769, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -7452.2014338237495, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 25, 'tm': 69.765715142391, 'gc_content': 0.6, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -6706.71143382373, 'success': True}
{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 24, 'tm': 68.83426028792064, 'gc_content': 0.625, 'hairpin_dg': -2383.1176867678478, 'homodimer_dg': -6391.038975986783, 'success': True}


