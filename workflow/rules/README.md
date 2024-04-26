
This workflow has rules / scripts for processing two types of data currently.
First I started with the more fine-grained and detailed base-by-base 
SMRF-seq calls. This data has the benefit of telling you exactly which bases
where converted but is more complex, the formating is more abstract and is in
general just harder to work with. Rules that were written to work with this data
are located in the `baseCallFiles.smk`

I soon switched over to just the bed files which only have the coordinates of
peaks on a per sample and plasmid basis. These are located in the 
`bedFiles.smk` rule file. 

To get everything organized for input into the snakemake pipeline for both
of these types of data took some organizing, namely merging peak calls that
were made for the same reads but for different conversion types. This took
place outside of the workflow and scripts to do this are located in the
`resources` dir. The reasoning for this was so I could produce nicer tsv files
with all wildcard variables I would need to it was easier to calculate exactly
what files were going into and out of the workflow which avoids having to
do things with snakemake checkpoints.