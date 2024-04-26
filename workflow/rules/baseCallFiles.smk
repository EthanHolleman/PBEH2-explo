
samples = pd.read_csv("../resources/samples/samples.grouped.tsv", sep="\t")

keep_rows = []

for i, each_row in samples.iterrows():
    with open(each_row.merge_path) as handle:
        n_lines = handle.readlines()
        if len(n_lines) > 2:
            keep_rows.append(i)

samples = samples.iloc[keep_rows]
groups = set(list(samples["group_id"]))


rule all:
    input:
        expand(
            "output/plots/peak.count.plots.pdf",
            group=groups,
        ),


rule make_composite_calls:
    conda:
        "envs/py.yml"
    input:
        lambda w: samples.loc[samples.group_id == int(w.group)].merge_path.to_list(),
    output:
        out="output/sampleCalls/composite.treat.type.group.{group}.tsv",
    script:
        "scripts/merge_conversions.py"


rule add_sample_treatments:
    conda:
        "envs/py.yml"
    resources:
        mem_mb=6000,
    input:
        sample_df="../resources/sampleTable.tsv",
        comp_df="output/sampleCalls/composite.treat.type.group.{group}.tsv",
    output:
        out="output/sampleCalls/composite.treat.group.{group}.tsv",
    script:
        "scripts/merge_sample_table.py"


rule add_call_type:
    conda:
        "envs/py.yml"
    input:
        tsv="output/sampleCalls/composite.treat.group.{group}.tsv",
    output:
        out="output/sampleCallsType/composite.treat.call.group.{group}.tsv",
    params:
        sample_table="../resources/samples/samples.grouped.tsv",
        group = lambda w: w.group
    script:
        "scripts/add_call_type.py"


rule prune_tsv_for_cat:
    conda:
        'envs/py.yml'
    input:
        tsv='output/sampleCallsType/composite.treat.call.group.{group}.tsv'
    output:
        out='output/sampleCallsTypePrune/composite.treat.call.prune.group.{group}.tsv'
    script:'scripts/prune_for_concat.py'
    


rule concat_peaks:
    conda:
        "envs/py.yml"
    resources:
        mem_mb=50000,
    input:
        tsvs=expand(
            "output/sampleCallsTypePrune/composite.treat.call.prune.group.{group}.tsv",
            group=groups,
        ),
    output:
        out="output/concatPeaks/peak.merge.tsv",
    script:
        "scripts/concat_tsvs.py"


rule plot_peak_counts:
    conda:
        "envs/R.yml"
    resources:
        mem_mb=50000,
    input:
        tsv="output/concatPeaks/peak.merge.tsv",
    output:
        out="output/plots/peak.count.plots.pdf",
    script:
        "scripts/plot.peak.counts.R"
