import pandas as pd


ID_VARS = [
    'sample_id',
    'start',
    'Experiment',
    'Experiment_name',
    'end',
    'call_type',
    'strand',
    'Plasmid',
    'Total plasmids in sample'
]

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


def main():

    df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    melted_df = df.melt(
        id_vars=ID_VARS, value_vars=VALUES, var_name='treatment',
        value_name='treatment_value'
    )
    melted_df.to_csv(snakemake.output['out'], sep='\t', index=False)


if __name__ == '__main__':
    main()
