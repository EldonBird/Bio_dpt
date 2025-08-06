import pandas

import pcr_lib as pl
import pandas as pd

data = {
    "snp_id": ["rs1001", "rs1002", "rs1003", "rs1004"],
    "sequence": ["ATCGA", "GGTAA", "CTGAA", "TACCA"],
    "allele": ["A", "G", "C", "T"],
    "position": [12345, 12367, 12400, 12450]
}

df = pandas.DataFrame(data)

print(df)


temp = pl.df_to_listprimers(df)

something = pl.generate_allele_specific_primers(temp, 1, 10)

for i in temp:
    print(i.sequence)

