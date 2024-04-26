# like the base by base calls, bed files also do not represent unique
# peaks as supplied because the same peak will be present for different
# types of C conversion types for both strands. For most cases I just care
# about having 1 file per strand that defines where the peaks are for a given
# call parameters. So here in this file I merge peaks in bed files by
# sample id, plasmid and strand. 

# When merging peaks first look by unique read name, peaks that are from the
# same read and have the same start and end should be compressed to one peak

import pandas as pd 

TEMP_HEADERS = [
    'description', 'start', 'end', 'read', 'score', 'strand', 'path'
]
GROUPS = ['sample_id', 'plasmid', 'strand', 'call_type']


def main():

    sample_table = pd.read_csv(
        'sampleTableBed.tsv',
        sep='\t'
    )

    group_list, names = make_group_list(sample_table)
    prune_samples = []
    for each_group, each_name in zip(group_list, names):
        peak_group = read_all_peaks(each_group)
        prune = prune_peaks(peak_group)
        prune_samples.append(
            write_prune_file(prune, each_name)
        )

    # Made the merged and pruned files now write a csv file telling me where
    # they are and how they were grouped for snakemake in case I want to use
    # that instead of wildcard approach

    pruned_samples_df = pd.DataFrame(prune_samples)
    pruned_samples_df.to_csv('samples/samples.bed.grouped.tsv', sep='\t', index=False)


def make_group_list(sample_table):

    grouped = sample_table.groupby(GROUPS)
    group_names = []
    group_list = []

    for name, group_df in grouped:
        group_list.append(group_df)
        group_names.append(name)
        

    return group_list, group_names


def read_all_peaks(group_df):
    bed_list = []
    for each_bed in group_df.filepath:
        bed_list.append(pd.read_csv(each_bed, sep='\t', header=None))
        bed_list[-1].columns = TEMP_HEADERS
    return pd.concat(bed_list)


def prune_peaks(peak_df):

    return peak_df.drop_duplicates(subset=['read', 'start', 'end'])


def write_prune_file(df, group_name):
   
    filepath = get_group_filepath(group_name)
    df.to_csv(filepath, sep='\t', index=False)
    
    sample_dict = dict(zip(GROUPS, group_name))
    sample_dict['filepath'] = filepath

    return sample_dict

def get_group_filepath(group_name):

    return f'mergedBed/PBEH2-{group_name[0]}_Plasmid_{group_name[1]}_Strand_{group_name[2]}_call_{group_name[3]}.bed.tsv'



if __name__ == '__main__':
    main()