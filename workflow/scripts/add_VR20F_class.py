# There are three classes of VR20F plasmid, high, med or low folding
# coeffiecients. This script just adds an additional column with this
# classification to make plotting and comparisons in R easier

import pandas as pd
import re

def main():

    df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    df['class'] = df.apply(
        lambda row: extract_plasmid_class(row),
        axis=1
    )
    df['version'] = df.apply(
        lambda row: extract_plasmid_class(row, group=1),
        axis=1
    )
    df.to_csv(snakemake.output['out'], sep='\t', index=False)


def extract_plasmid_class(row, group=0):

    match = re.search(r'VR20F_\d(.+)(\d)', row.plasmid)
    return match.groups()[group]


if __name__ == '__main__':
    main()




