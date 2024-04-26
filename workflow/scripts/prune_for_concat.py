import pandas as pd

# Get only the essential columns for plotting to avoid creating huge
# tsv file and slowing things down

VALUES = [
    'Linear',
    'Nick',
    'Txn',
    'RnaseH',
    'Co-SSB',
    'Post-SSB',
    'Co-RPA',
    'Linearization post txn',
    'Deprotienize post Co-RPA'
]

TARGET_COLS = [
    'sample', 'sample_id', 'start', 'end', 'Experiment', 'Experiment_name',
    'strand', 'call_type', 'PEAK_x'
]

KEEP = VALUES + TARGET_COLS

def main():


    df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    df_keep = df[KEEP]
    df_keep.to_csv(snakemake.output['out'], sep='\t', index=False)


if __name__ == '__main__':
    main()