
import pandas as pd


def main():
    
    sample_labels = pd.read_csv(snakemake.input['sample_labels'], sep='\t')
    peak_freq = pd.read_csv(snakemake.input['tsv'], sep='\t')

    sample_labels = get_int_sample_id(sample_labels)
    
    merge = peak_freq.merge(sample_labels, on='sample_id')
    
    merge.to_csv(snakemake.output['out'], sep='\t', index=False)



def get_int_sample_id(sample_labels):
    
    def extract_sample_id(row):
        
        return int(row.Sample_ID.split('-')[-1])
    
    sample_labels['sample_id'] = sample_labels.apply(
        lambda row: extract_sample_id(row),
        axis=1
    )
    
    return sample_labels


if __name__ == '__main__':
    main()