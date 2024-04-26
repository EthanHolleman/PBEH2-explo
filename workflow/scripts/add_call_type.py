import pandas as pd


def main():

    #df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    #sample_table = pd.read_csv(snakemake.params['sample_table'], sep='\t')
    #sample_table_call_type = sample_table[['sample_id', 'call_type']]

    df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    st = pd.read_csv(snakemake.params['sample_table'], sep='\t')
    group = int(snakemake.params['group'])
    
    call_type = (set(st.loc[st.group_id==group].call_type))

    assert len(call_type) == 1  # should only be one call type since its by group

    df['call_type'] = call_type.pop()

    df.to_csv(snakemake.output['out'], sep='\t', index=False)

if __name__ == '__main__':
    main()