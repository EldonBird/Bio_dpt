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
            alleles = allele_str.split("/") if allele_str else []

            # Ensure there are at least two alleles to work with
            if len(alleles) < 2:
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
    # snp_df = Fetch_SNP_Data(["rs1799971"], 30)# just here for testing.  , "rs599839"
  

    # sudo data to speed up testing
    snp_df = [{'snpID': 'rs1799971', 'allele': 'A', 'sequence': 'TGTGTTTGCACAGAAGAGTGCCCAGTGAAGAGACCTACTCCTTGGATCGCTTTGCGCAAAATCCACCCCTTTTCCCTCCTCCCTCCCTTCCAGCCTCCGAATCCCGCATGGCCCACGCTCCCCTCCTGCAGCGGTGCGGGGCAGGTGATGAGCCTCTGTGAACTACTAAGGTGGGAGGGGGCTATACGCAGAGGAGAATGTCAGATGCTCAGCTCGGTCCCCTCCGCCTGACGCTCCTCTCTGTCTCAGCCAGGACTGGTTTCTGTAAGAAACAGCAGGAGCTGTGGCAGCGGCGAAAGGAAGCGGCTGAGGCGCTTGGAACCCGAAAAGTCTCGGTGCTCCTGGCTACCTCGCACAGCGGTGCCCGCCCGGCCGTCAGTACCATGGACAGCAGCGCTGCCCCCACGAACGCCAGCAATTGCACTGATGCCTTGGCGTACTCAAGTTGCTCCCCAGCACCCAGCCCCGGTTCCTGGGTCAACTTGTCCCACTTAGATGGCAACCTGTCCGACCCATGCGGTCCGAACCGCACCGACCTGGGCGGGAGAGACAGCCTGTGCCCTCCGACCGGCAGTCCCTCCATGATCACGGCCATCACGATCATGGCCCTCTACTCCATCGTGTGCGTGGTGGGGCTCTTCGGAAACTTCCTGGTCATGTATGTGATTGTCAGGTAAGGAAAGCGCCAGGGCTCCGAGCGGAGGGTTCAGCGGCTTAAGGGGGTACAAAGAGACACCTAACTCCCAAGGCTCAATGTTGGGCGGGAGGATGAAAGAGGGGAGGTAAACTGGGGGGACTCTGGAGGAGACCACGGACAGTGATTGTTATTTCTATGAGAAAACCTACTTTTCTGTTTTTTCTTCAACTGATAAAGAAAGAATTCAAAATTTCAGGAGCAGAGAAGTTGCTTTGGTAAAAGCTACAAATGTCTAGGGGTGGGGGGCGGAGGGAAGCTATAGCATAGACTTGGAGCGCTTCCTTATACTGAGCAAAGAGGGCTC', 'position': 500}, {'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TGTGTTTGCACAGAAGAGTGCCCAGTGAAGAGACCTACTCCTTGGATCGCTTTGCGCAAAATCCACCCCTTTTCCCTCCTCCCTCCCTTCCAGCCTCCGAATCCCGCATGGCCCACGCTCCCCTCCTGCAGCGGTGCGGGGCAGGTGATGAGCCTCTGTGAACTACTAAGGTGGGAGGGGGCTATACGCAGAGGAGAATGTCAGATGCTCAGCTCGGTCCCCTCCGCCTGACGCTCCTCTCTGTCTCAGCCAGGACTGGTTTCTGTAAGAAACAGCAGGAGCTGTGGCAGCGGCGAAAGGAAGCGGCTGAGGCGCTTGGAACCCGAAAAGTCTCGGTGCTCCTGGCTACCTCGCACAGCGGTGCCCGCCCGGCCGTCAGTACCATGGACAGCAGCGCTGCCCCCACGAACGCCAGCAATTGCACTGATGCCTTGGCGTACTCAAGTTGCTCCCCAGCACCCAGCCCCGGTTCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCACCGACCTGGGCGGGAGAGACAGCCTGTGCCCTCCGACCGGCAGTCCCTCCATGATCACGGCCATCACGATCATGGCCCTCTACTCCATCGTGTGCGTGGTGGGGCTCTTCGGAAACTTCCTGGTCATGTATGTGATTGTCAGGTAAGGAAAGCGCCAGGGCTCCGAGCGGAGGGTTCAGCGGCTTAAGGGGGTACAAAGAGACACCTAACTCCCAAGGCTCAATGTTGGGCGGGAGGATGAAAGAGGGGAGGTAAACTGGGGGGACTCTGGAGGAGACCACGGACAGTGATTGTTATTTCTATGAGAAAACCTACTTTTCTGTTTTTTCTTCAACTGATAAAGAAAGAATTCAAAATTTCAGGAGCAGAGAAGTTGCTTTGGTAAAAGCTACAAATGTCTAGGGGTGGGGGGCGGAGGGAAGCTATAGCATAGACTTGGAGCGCTTCCTTATACTGAGCAAAGAGGGCTC', 'position': 500}]


    primers = generate_allele_specific_primers(snp_df, 24, 26)
    # print(primers)
    print(f"number of snps {len(primers)}")
    print(f"number of directions {len(primers[0])}")
    print(f"number of primers in a direction {len(primers[0][0])}")
    print(len(snp_df))
    # for snp in primers:



if(__name__ == "__main__"):
    Main()


# [[[{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TTAGATGACA', 'direction': 'forward', 'length': 10, 'tm': 18.129666660821954, 'gc_content': 0.3, 'hairpin_dg': 0.0, 'homodimer_dg': -72.98389370605668, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TAGATGACA', 'direction': 'forward', 'length': 9, 'tm': 12.174837327327339, 'gc_content': 0.3333333333333333, 'hairpin_dg': 0.0, 'homodimer_dg': -72.98389370605668, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'AGATGACA', 'direction': 'forward', 'length': 8, 'tm': 7.16645091498475, 'gc_content': 0.375, 'hairpin_dg': 0.0, 'homodimer_dg': -72.98389370605668, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GATGACA', 'direction': 'forward', 'length': 7, 'tm': -1.3545307221136227, 'gc_content': 0.42857142857142855, 'hairpin_dg': 0.0, 'homodimer_dg': -72.98389370605668, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'ATGACA', 'direction': 'forward', 'length': 6, 'tm': -17.656890150261773, 'gc_content': 0.3333333333333333, 'hairpin_dg': 0.0, 'homodimer_dg': -72.98389370605668, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TGACA', 'direction': 'forward', 'length': 5, 'tm': -32.72779405071134, 'gc_content': 0.4, 'hairpin_dg': 0.0, 'homodimer_dg': -90.36415040613073, 'success': True}], 
   
#   [{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGACAGATT', 'direction': 'reverse', 'length': 10, 'tm': 27.203338243722783, 'gc_content': 0.5, 'hairpin_dg': 0.0, 'homodimer_dg': -73.96051300974477, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGACAGATT', 'direction': 'reverse', 'length': 9, 'tm': 16.740737528431566, 'gc_content': 0.4444444444444444, 'hairpin_dg': 0.0, 'homodimer_dg': -73.96051300974477, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GACAGATT', 'direction': 'reverse', 'length': 8, 'tm': 6.159145398662417, 'gc_content': 0.375, 'hairpin_dg': 0.0, 'homodimer_dg': -73.96051300974477, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'ACAGATT', 'direction': 'reverse', 'length': 7, 'tm': -6.651097932493769, 'gc_content': 0.2857142857142857, 'hairpin_dg': 0.0, 'homodimer_dg': -73.96051300974477, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CAGATT', 'direction': 'reverse', 'length': 6, 'tm': -19.54844111334802, 'gc_content': 0.3333333333333333, 'hairpin_dg': 0.0, 'homodimer_dg': -114.00720323817018, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'AGATT', 'direction': 'reverse', 'length': 5, 'tm': -45.52810475874327, 'gc_content': 0.2, 'hairpin_dg': 0.0, 'homodimer_dg': -130.5585090884421, 'success': True}]]]


# [[[{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GGTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 26, 'tm': 62.64091186246884, 'gc_content': 0.46153846153846156, 'hairpin_dg': 40.246484286148586, 'homodimer_dg': 11.618441795495414, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'GTCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 25, 'tm': 60.741963776567445, 'gc_content': 0.44, 'hairpin_dg': 39.085883205265475, 'homodimer_dg': -11.318504829494259, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCAACTTGTCCCACTTAGATGACA', 'direction': 'forward', 'length': 24, 'tm': 59.65329123118744, 'gc_content': 0.4166666666666667, 'hairpin_dg': 36.03872822148651, 'homodimer_dg': -32.43413848254389, 'success': True}], 
  
#   [{'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TTCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 26, 'tm': 69.86363516091336, 'gc_content': 0.5769230769230769, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 27.76279186528285, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'TCGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 25, 'tm': 69.765715142391, 'gc_content': 0.6, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 25.839909279650783, 'success': True}, 
#    {'snpID': 'rs1799971', 'allele': 'A', 'primer_sequence': 'CGGACCGCATGGGTCGGACAGATT', 'direction': 'reverse', 'length': 24, 'tm': 68.83426028792064, 'gc_content': 0.625, 'hairpin_dg': 61.78877378212809, 'homodimer_dg': 24.30006522953198, 'success': True}]]]