# Add columns from the sample table data like treatments and experiment names
# to the grouped read counts table for easier plotting by experiment name

import pandas as pd


def main():
    
    
    sample_treatments = pd.read_csv(snakemake.input['sample_treatments'], sep='\t')
    sample_treatments = add_integer_sample_id(sample_treatments)
    
    sample_counts = pd.read_csv(snakemake.input['sample_groups'], sep='\t')
    
    merge = sample_counts.merge(sample_treatments, on='sample_id')
    
    merge.to_csv(snakemake.output['out'], sep='\t', index=False)


def add_integer_sample_id(sample_treatments):
    
    def extract_int_sample_id(row):
        return int(row['Sample_ID'].split('-')[-1])
    
    sample_treatments['sample_id'] = sample_treatments.apply(
        lambda row: extract_int_sample_id(row),
        axis=1
    )
    return sample_treatments


if __name__ == '__main__':
    main()