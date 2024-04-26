import pandas as pd

df = pd.read_csv('/home/ethollem/workflows/PBEH2-explo/workflow/output/bedData/mergedExperiments/Nt.BspqI.nicked.samples.named.tsv', sep='\t')

read_count = df.loc[
    (df['call_type'] == 'low') & 
    (df['strand'] == 'Pos') &
    (df['plasmid'] == 'PFC9_NTBSPQI_11') &
    (df['treatment_phrase'] == 'Linear Txn-')
    
]

read_count.to_csv('test.count.tsv', sep='\t')