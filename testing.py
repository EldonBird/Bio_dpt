from Bio import Seq
import numpy as np
import prc_lib as plib
# snp = 'ATGCAATTGGCCAAATTTGGGCCCAAAATTTTGGGGCCCCAAAAATTTTTGGGGGCCCCCAAAAAATTTTTTGGGGGGCCCCCC'
# max_len = 30
# min_len = 18
# for length in range(max_len-min_len-1):
#     trimmed = snp[length]
#     print(trimmed)
#     print(length)


snp_data = [{'snpID': 'rs1799971', 'allele': 'A', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCAACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, {'snpID': 'rs1799971', 'allele': 'G', 'sequence': 'TCCTGGGTCAACTTGTCCCACTTAGATGGCGACCTGTCCGACCCATGCGGTCCGAACCGCA', 'position': 30}, {'snpID': 'rs599839', 'allele': 'G', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCGACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, {'snpID': 'rs599839', 'allele': 'A', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCAACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, {'snpID': 'rs599839', 'allele': 'C', 'sequence': 
'AAAAAAAGAGAAAGAAATAGGAGCAGGATCCACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}, {'snpID': 'rs599839', 'allele': 'T', 'sequence': 'AAAAAAAGAGAAAGAAATAGGAGCAGGATCTACTTCCAGATATACAGAGAATATAAAAATA', 'position': 30}]

# snp_short = snp[0:10]
# print(f"Normal string: {snp_short}")
# seq = Seq.Seq(snp_short)
# print(f"Now a sequence in BIO Python: {seq}")
# reverse = seq.reverse_complement()
# print(f"Now a reverse sequence in Bio Python: {reverse}")
# rev_string = str(reverse)
# print(f"Now a string again: {rev_string}")
# reverse

# np_frame = np.array(snp_data)

# for i in np_frame:
#     print(i['sequence'])
plib.generate_allele_specific_primers(snp_data, 24, 26)