import pandas as pd
import re

READ_COUNT_TSV = [
    'PBEH2-revio-read-counts/PBEH2_240409_L1_t45_w15_l50_Neg_GH_stat.tsv',
    'PBEH2-revio-read-counts/PBEH2_240409_L1_t45_w15_l50_Pos_CH_stat.tsv',
    'PBEH2-revio-read-counts/PBEH2_240409_L1_t55_w20_l100_Neg_GH_stat.tsv',
    'PBEH2-revio-read-counts/PBEH2_240409_L1_t55_w20_l100_Pos_CH_stat.tsv'
]

READ_COUNT_DFS = [pd.read_csv(tsv, sep='\t') for tsv in READ_COUNT_TSV]


def main():

    sample_table = pd.read_csv('samples/samples.bed.grouped.tsv', sep='\t')

    # Pull out the interger sample id so can look up by sample easier and
    # the plasmid name
    read_count_dfs = [
        add_plasmid_to_read_count_df(add_sample_id_to_read_count_df(df)) 
        for df in READ_COUNT_DFS
        ]
    
    # Add the read count to each sample
    sample_table['read_count'] = sample_table.apply(
        lambda row: assign_read_count(row),
        axis=1
    )

    sample_table.to_csv(
        'samples/samples.bed.grouped.count.tsv', sep='\t', index=False
        )

    

def pick_read_count_df(strand, call_type):

    if strand == 'Pos':

        if call_type == 'high':
            return READ_COUNT_DFS[3]
        else:
            return READ_COUNT_DFS[1]
    
    else:
        if call_type == 'high':
            return READ_COUNT_DFS[2]
        else:
            return READ_COUNT_DFS[0]


def add_sample_id_to_read_count_df(read_count_df):

    def extract_sample_id(sample_name):
        match = re.search(r'BCBC(\d+)', sample_name)
        sample_id = match.groups()[0]
        if sample_id[0] == '0':
            sample_id = int(sample_id[1])
        else:
            sample_id = int(sample_id)
        return sample_id
    
    read_count_df['sample_id'] = read_count_df.apply(
        lambda row: extract_sample_id(row['#sample']),
        axis=1
    )

    return read_count_df


def add_plasmid_to_read_count_df(read_count_df):

    def extract_plasmid_name(sample_name):
        match = re.search(r'PLASMID(.+)_DESC', sample_name)
        return match.groups()[0]
    

    read_count_df['plasmid'] = read_count_df.apply(
        lambda row: extract_plasmid_name(row['#sample']),
        axis=1
    )

    return read_count_df




def assign_read_count(row):
    
    read_count_df = pick_read_count_df(row.strand, row.call_type)

    read_count = read_count_df.loc[
        (read_count_df['sample_id'] == row.sample_id) & 
        (read_count_df['strand'] == row.strand) &
         (read_count_df['plasmid'] == row.plasmid)
    ]

    # special case for PFC8_TACINIT series since they had to be "fixed" but
    # the plasmid name was not changed in a way that would adhere to the style
    # of the rest of the file formating

    samples = list(read_count['#sample'])
    for each_sample in samples:
        if 'PFC8_TACINIT' in each_sample:
            # Use only the "fixed" read counts. Add a column 
            read_count['fixed'] = read_count.apply(
                lambda row: row['#sample'].split('_')[-1],
                axis=1
            )
            read_count = read_count.loc[read_count.fixed=='FIXED']
            break
    
    assert len(read_count) == 1, f'Should only have 1 read count, you had {len(read_count)}'

    return list(read_count['total_read'])[0]


if __name__ == '__main__':
    main()