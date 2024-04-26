import pandas as pd


def main():

    tsv_files = snakemake.input['tsvs']
    dfs = [pd.read_csv(t, sep='\t') for t in tsv_files]
    concat_df = pd.concat(dfs)
    concat_df.to_csv(snakemake.output['out'], sep='\t', index=False)


if __name__ == '__main__':
    main()
