
# When getting sample frequency data it comes labeled with a table of what
# treatments where applied to each sample. Each treatment is a column and if
# it was applied the value for each row will be + or - or a specific treatment
# type like an enzyme. The purpose of this script is to turn those treatment
# tables into a short phrase that sums up what that sample treatment represents
# and create a new column with this phrase so it can be used for plotting

import pandas as pd


def main():
    
    df = pd.read_csv(snakemake.input['tsv'], sep='\t')
    
    # Figure out if any of the samples have a co or post ssb treatment
    # if at least one does then treat this experiment as an ssb experiment
    # this is done so we dont have to add no ssb to experiments where nothing
    # was treated with SSB. In the future do the same thing with RPA
    
    ssb_experiment = is_ssb_experiment(df)
    
    df['treatment_phrase'] = df.apply(
        lambda row: determine_treatment_phrase(row, ssb_experiment),
        axis=1
    )
    df.to_csv(snakemake.output['out'], sep='\t', index=False)


def determine_treatment_phrase(row, ssb_experiment=False):
    
    phrase = []
    
    if row.Nick != '-':
        phrase.append('Nicked')
    
    if row.Linear != '-':
        phrase.append('Linear')
    else:
        phrase.append('SC')
        
    if ssb_experiment:
        if row['Co-SSB'] == '+':
            phrase.append('Co-SSB')
        elif row['Post-SSB'] == '+':
            phrase.append('Post-SSB')
        else:
            phrase.append('No SSB')
        
        if row['Linearization post txn'] == '+':
            phrase.append('Cut post txn')
    
    if row.Txn == '+':  # sample was transcribed
        if row.RnaseH == '+':
            phrase.append('RnaseH')
        else:
            phrase.append('Txn+')
    else:
        phrase.append('Txn-')
    
    return ' '.join(phrase)


def is_ssb_experiment(df):
    
    if '+' in set(df['Co-SSB']) or '+' in set(df['Post-SSB']):
        return True
    else:
        return False


if __name__ == '__main__':
    main()
        
        