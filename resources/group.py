import pandas as pd
import re



def main():

    samples = pd.read_csv('/home/ethollem/workflows/rloop_fold/workflow/samples/samples.tsv', sep='\t')
    # Add plasmid name

    samples['plasmid'] = samples.apply(
        lambda row: extract_plasmid_name(row),
        axis=1
    )

    grouped = samples.groupby(['call_type', 'sample_id', 'strand', 'file_type', 'plasmid'])

    dfs = [grouped.get_group(group) for group in grouped.groups]

    for i, each_df in enumerate(dfs):
        each_df['group_id'] = i

    with_groups = pd.concat(dfs)

    with_groups.to_csv('samples.grouped.tsv', sep='\t', index=False)



def extract_plasmid_name(row):
    """Get the name of the plasmid from the sample description

    Args:
        row (Series): Peak dataframe row

    Returns:
        str: Name of the plasmid peak is called on
    """

    plasmid = re.search(r'(.+)_DESC', row['description']).groups()[0]
    return plasmid



if __name__ == '__main__':
    main()
