import pandas as pd

# All frames should have the same peaks called in them so use the first
# frame that we read to get the peak names


# We really just care about Cs that are called as peaks and those that
# are not. 7 or 8 means called in a peak so for each read get all of the
# different call types from the different files and make a "consensus" call
# which marks 7 or 8 if any of the calls have 7 or 8 in that position

def make_composite_call(call_types):


    def extract_read_seq(peak_row):

        colnames = peak_row.index[14:-1]
        seq = [col.split('.')[0] for col in colnames]
        return ''.join(seq)


    def get_calls_from_frame(peak_row):

        # Sequence starts at column 14

        return list(peak_row[14:-1])
    
    def get_id_info_from_frame(peak_row):

        info = peak_row[:13]
        names = peak_row.index[:13]
        return dict(zip(names, info))
    
    # get all of the calls as a list
    calls = [get_calls_from_frame(peak_row) for _, peak_row in call_types.iterrows()]
    info = get_id_info_from_frame(call_types.iloc[0])
    seq = extract_read_seq(call_types.iloc[0])

    # use -1 as place holder values
    consensus_call = [-1] * len(calls[0])
  
    for each_call_list in calls:
        for i, each_call in enumerate(each_call_list):
            # if value has not been updated (-1) then update it
            if consensus_call[i] == -1:
                consensus_call[i] = each_call
            
            # Only overwrite existing values if the call is for a converted
            # base
            if each_call == 8 or each_call == 9:
                consensus_call[i] = each_call
    
    consensus_call = ''.join([str(i) for i in consensus_call])

    info['composite_call'] = consensus_call
    info['sequence'] = seq

    return info


def main():

    # files from the same sample, strand and call type
    sample_files = snakemake.input
    print(sample_files)
    #sample_files = '../resources/mergedFiles/PBEH2_bismark_bt2.bam_genePBEH2_BCBC91_PLASMIDT7INIT_VR27_DESCSUPERCOIL_NT_Neg_15_0.45_GC.PEAK.genome.merge, ../resources/mergedFiles/PBEH2_bismark_bt2.bam_genePBEH2_BCBC91_PLASMIDT7INIT_VR27_DESCSUPERCOIL_NT_Neg_15_0.45_CG.PEAK.genome.merge, ../resources/mergedFiles/PBEH2_bismark_bt2.bam_genePBEH2_BCBC91_PLASMIDT7INIT_VR27_DESCSUPERCOIL_NT_Neg_15_0.45_CH.PEAK.genome.merge, ../resources/mergedFiles/PBEH2_bismark_bt2.bam_genePBEH2_BCBC91_PLASMIDT7INIT_VR27_DESCSUPERCOIL_NT_Neg_15_0.45_GH.PEAK.genome.merge'
    if isinstance(sample_files, str):  # single file case
        sample_files = [sample_files]


    frames = [pd.read_csv(sample, sep='\t') for sample in sample_files]

    if len(frames) > 1:
        cat_frames = pd.concat(frames)
    else:
        cat_frames = frames[0]

    peaks = set(list(cat_frames.PEAK_x))

    composite_calls = []

    for each_peak in peaks:

        call_types = cat_frames.loc[cat_frames.PEAK_x == each_peak]
        composite_calls.append(
            make_composite_call(call_types)
        )
        assert isinstance(composite_calls[-1]['composite_call'], str)
    
    assert len(composite_calls) == len(peaks)
    
    c_calls_df = pd.DataFrame(composite_calls)
    c_calls_df.to_csv(snakemake.output['out'], sep='\t', index=False)



if __name__ == '__main__':
    main()
    




