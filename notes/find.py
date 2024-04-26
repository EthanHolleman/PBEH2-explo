import pandas as pd

d = '/home/ethollem/workflows/PBEH2-explo/resources/samples/samples.bed.grouped.count.tsv'
df = pd.read_csv(d, sep='\t')
#  plasmid=PFC9_NTBSPQI_08, strand=Pos, call_type=high
s = df.loc[(df.sample_id == 28) &
                    (df.plasmid == 'PFC9_NTBSPQI_08') &
                    (df.call_type == 'high') &
                   (df.strand == 'Pos')
                    ]

print(s)