# Iterate through SMRF-seq bed files and create a dataframe with sample
# information for each file (plasmid, sample id, etc) as well as the orignal
# bed file path so can be more easily accessed in snakemake


import re
from pathlib import Path
import pandas as pd


LOW_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t45_w15_l50/PEAKS_GENOME'
HIGH_DIR = '/group/flchedingrp/mitochi/Work/Project/Ethan/240325/5_Analysis/2_deleted/67_72_PFC8TACINIT_240406/PBEH2_240409_L1_t55_w20_l100/PEAKS_GENOME'



def main():

    bed_files_high = get_bed_files_from_dir(HIGH_DIR)
    bed_files_low = get_bed_files_from_dir(LOW_DIR)

    high_calls = [
        create_sample_dict(each_bed, 'high') for each_bed in bed_files_high
    ]
    low_calls = [
        create_sample_dict(each_bed, 'low') for each_bed in bed_files_low
    ]

    assert len(high_calls) == len(low_calls), f'{len(high_calls)}-{len(low_calls)}'
    assert len(high_calls) == 5118

    all_calls = high_calls + low_calls

    # add a unique number for each bed file and correct sample id

    for i, each_bed in enumerate(all_calls):
        each_bed['bed_id'] = i
        each_bed['sample_id'] = convert_sample_id_to_int(each_bed['sample_id'])
    
    bed_df = pd.DataFrame(all_calls)
    bed_df = bed_df.sort_values(['sample_id', 'call_type', 'plasmid', 'strand'])
    bed_df.to_csv('sampleTableBed.tsv', sep='\t', index=False)


def convert_sample_id_to_int(sample_id_str):

    if sample_id_str[0] == '0':
        sample_id = sample_id_str[1:]
    else:
        sample_id = sample_id_str

    sample_id = int(sample_id)

    return sample_id


def create_sample_dict(filepath, call_type):

    def get_sample_info_from_filename():

        filename = filepath.name
        m = re.search(
            r'BCBC(\d+)_PLASMID(.+)_DESC(.+)_(Pos|Neg)_(.+)_(.+).PEAK.genome.bed',
            filename
        )
        return dict(zip(
            ['sample_id', 'plasmid', 'description', 'strand', 'window', 'convert_type'],
            m.groups())
        )

    sample_dict = get_sample_info_from_filename()
    sample_dict['filepath'] = str(filepath)
    sample_dict['call_type'] = call_type

    return sample_dict

def get_sample_info_from_filename(filename):

    m = re.search(
        r'BCBC(\d+)_PLASMID(.+)_DESC(.+)_(Pos|Neg)_(.+)_(.+).PEAK.genome.bed'
    )
    return dict(zip(
        ['sample_id', 'plasmid', 'description', 'strand', 'window', 'convert_type'],
        m.groups()
    ))


def get_bed_files_from_dir(dir):
    paths = []
    for each_path in Path(dir).iterdir():
        if each_path.suffix == '.bed':
            paths.append(each_path)
    return paths


if __name__ == '__main__':
    main()