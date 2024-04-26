import pandas as pd
import numpy as np


# bed files are seperated out by plasmid, sample and strand so now just need
# to compile the peaks and divide by read length but the next sad thing
# is that I need to know the plasmid length, probably would be good to determine
# that prior to this step but for now just create a list long enough to accomidate
# all the conversion locations


def main():

    peaks = snakemake.input['tsv']
    samples = snakemake.input['sample_table']
    #peaks = '../resources/mergedBed/PBEH2-33_Plasmid_VR20F_1LOW1_Strand_Neg_call_low.bed.tsv'
    #samples = '../resources/samples/samples.bed.grouped.count.tsv'
    
    peak_df = pd.read_csv(peaks, sep='\t')
    sample_df = pd.read_csv(samples, sep='\t')

    read_len = 3500
    total_reads = 2500 #get_total_reads(sample_df)

    freq_df = make_frequency_length_df(peak_df, total_reads, read_len)

    freq_df.to_csv(snakemake.output['out'], sep='\t')


def get_total_reads(df: pd.DataFrame):

    sample = df.loc[(df.sample_id == int(snakemake.wildcards.sample_id)) &
                    (df.strand == snakemake.wildcards.strand) & 
                    (df.call_type == snakemake.wildcards.call_type) &
                    (df.plasmid == snakemake.wildcards.plasmid)
                    ]
    assert len(sample) == 1, sample

    return list(sample['read_count'])[0]


def make_frequency_length_df(peaks_df: pd.DataFrame, total_reads: int, read_len: int):
    '''Calculate the number of peaks, their frequency and the mean peak length
    at each base pair of a read.

    Args:
        peaks_df (pd.DataFrame): Dataframe with peak info (bed file)
        total_reads (int): Total number of reads for this sample
        read_len (int): Expected mean length 

    Returns:
        _type_: pd.DataFrame
    '''
    counts = np.zeros(read_len)
    
    num_peaks = len(peaks_df)
    lengths =  [[] for i in range(read_len)]

    for i, each_peak in peaks_df.iterrows():
        counts[each_peak.start:each_peak.end] += 1
        
        peak_len = each_peak.end - each_peak.start
        
        for each_bp in range(each_peak.start, each_peak.end):
            lengths[each_bp].append(peak_len)
        
    
    count_list = []
    
    def get_length_stats():
        
        length_stats = []
        for each_base in lengths:
            #print('MEAN', np.mean(each_base))
            if len(each_base) > 1:
                mean_length = np.mean(each_base)
                std = np.std(each_base)
            elif len(each_base) == 1:
                mean_length = each_base[0]
                std = 0
            else:
                mean_length = 0
                std = 0
            
            length_stats.append(
                {
                    'mean_len': mean_length,
                    'std': std
                }
            )
        
        return length_stats


    length_stats = get_length_stats()
    print(lengths[0])
    
    for i, each_count in enumerate(counts):

        count_dict = {
            'position': i,
            'total_peaks': each_count,
            'peak_freq': round(each_count / total_reads, 3),
            'plasmid': snakemake.wildcards.plasmid,
            'sample_id': int(snakemake.wildcards.sample_id),
            'strand': snakemake.wildcards.strand,
            'call_type': snakemake.wildcards.call_type
        }
        
        count_dict.update(length_stats[i])

        count_list.append(count_dict)
    
    return pd.DataFrame(count_list)

if __name__ == '__main__':
    main()








