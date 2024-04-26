# Experiments of interest

BED_SAMPLES = pd.read_csv(
    "../resources/samples/samples.bed.grouped.count.tsv", sep="\t"
)

NTBSPQI_SAMPLES = BED_SAMPLES.loc[BED_SAMPLES.sample_id.isin(set(range(0, 18)))]
VR20F_SAMPLES = BED_SAMPLES.loc[BED_SAMPLES.sample_id.isin(set(range(30, 36)))]
NTBSPQI_SSB_SAMPLES = BED_SAMPLES.loc[BED_SAMPLES.sample_id.isin(set(range(18, 30)))]


rule calculate_base_conversion_freq:
    conda:
        "../envs/py.yml"
    input:
        tsv="../resources/mergedBed/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.bed.tsv",
        sample_table="../resources/samples/samples.bed.grouped.count.tsv",
    output:
        out="output/bedData/basepairConversionFreq/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.tsv",
    script:
        "../scripts/calc_convert_freq.py"


rule merge_sample_treatments:
    conda:
        "../envs/py.yml"
    input:
        tsv="output/bedData/basepairConversionFreq/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.tsv",
        sample_labels="../resources/sampleTable.tsv",
    output:
        out="output/bedData/basepairConversionFreqSampleLabels/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.label.tsv",
    script:
        "../scripts/merge_bed_samps.py"


rule merge_ntbpsqi_experiments:
    conda:
        "../envs/py.yml"
    input:
        tsvs=expand(
            "output/bedData/basepairConversionFreqSampleLabels/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.label.tsv",
            zip,
            sample_id=NTBSPQI_SAMPLES.sample_id,
            plasmid=NTBSPQI_SAMPLES.plasmid,
            strand=NTBSPQI_SAMPLES.strand,
            call_type=NTBSPQI_SAMPLES.call_type,
        ),
    output:
        out="output/bedData/mergedExperiments/Nt.BspqI.nicked.samples.tsv",
    script:
        "../scripts/concat_tsvs.py"


rule merge_VR20F_experiments:
    conda:
        "../envs/py.yml"
    input:
        tsvs=expand(
            "output/bedData/basepairConversionFreqSampleLabels/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.label.tsv",
            zip,
            sample_id=VR20F_SAMPLES.sample_id,
            plasmid=VR20F_SAMPLES.plasmid,
            strand=VR20F_SAMPLES.strand,
            call_type=VR20F_SAMPLES.call_type,
        ),
    output:
        out="output/bedData/mergedExperiments/VR20F.samples.tsv",
    script:
        "../scripts/concat_tsvs.py"


rule merge_NT_SSB_experiments:
    conda:
        "../envs/py.yml"
    input:
        tsvs=expand(
            "output/bedData/basepairConversionFreqSampleLabels/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.label.tsv",
            zip,
            sample_id=NTBSPQI_SSB_SAMPLES.sample_id,
            plasmid=NTBSPQI_SSB_SAMPLES.plasmid,
            strand=NTBSPQI_SSB_SAMPLES.strand,
            call_type=NTBSPQI_SSB_SAMPLES.call_type,
        ),
    output:
        out="output/bedData/mergedExperiments/Nt.BspQI.SSB.samples.tsv",
    script:
        "../scripts/concat_tsvs.py"


rule add_ntbspqI_experiment_sample_labels:
    conda:
        "../envs/py.yml"
    input:
        tsv="output/bedData/mergedExperiments/Nt.BspqI.nicked.samples.tsv",
    output:
        out="output/bedData/mergedExperiments/Nt.BspqI.nicked.samples.named.tsv",
    script:
        "../scripts/add_sample_type.py"


rule add_vr20F_experiment_sample_labels:
    conda:
        "../envs/py.yml"
    input:
        tsv="output/bedData/mergedExperiments/VR20F.samples.tsv",
    output:
        out="output/bedData/mergedExperiments/VR20F.samples.named.tsv",
    script:
        "../scripts/add_sample_type.py"


rule add_nt_ssb_experiment_sample_labels:
    conda:
        "../envs/py.yml"
    input:
        tsv="output/bedData/mergedExperiments/Nt.BspQI.SSB.samples.tsv",
    output:
        out="output/bedData/mergedExperiments/Nt.BspQI.SSB.samples.named.tsv",
    script:
        "../scripts/add_sample_type.py"


rule add_vr20_plasmid_class:
    conda:
        "../envs/py.yml"
    input:
        tsv="output/bedData/mergedExperiments/VR20F.samples.named.tsv",
    output:
        out="output/bedData/mergedExperiments/VR20F.samples.named.class.tsv",
    script:
        "../scripts/add_VR20F_class.py"


rule plot_ntpspqi_experiments:
    conda:
        "../envs/R.yml"
    input:
        tsv="output/bedData/mergedExperiments/Nt.BspqI.nicked.samples.named.tsv",
    output:
        out="output/plots/Nt.BspQI.experiment.freq.plots.pdf",
    script:
        "../scripts/plots/plot.ntbpsqi.R"


rule plot_vr20F_experiments:
    conda:
        "../envs/R.yml"
    input:
        tsv="output/bedData/mergedExperiments/VR20F.samples.named.class.tsv",
    output:
        out="output/plots/VR20F.experiment.freq.plots.pdf",
    script:
        "../scripts/plots/plot.vr20F.R"


rule plot_NT_SSB_experiments:
    conda:
        "../envs/R.yml"
    input:
        tsv="output/bedData/mergedExperiments/Nt.BspQI.SSB.samples.named.tsv",
    output:
        out="output/plots/Nt.BspQI.SSB.experiment.freq.plots.pdf",
    script:
        "../scripts/plots/plot.nt.ssb.R"


rule add_treatments_to_sample_groups_table:
    conda:
        "../envs/py.yml"
    input:
        sample_treatments="../resources/sampleTable.tsv",
        sample_groups="../resources/samples/samples.bed.grouped.count.tsv",
    output:
        out="../resources/samples/samples.bed.grouped.count.treatment.tsv",
    script:
        "../scripts/add_sample_info_to_read_counts.py"


rule add_treatments_phrase_groups_table:
    conda:
        "../envs/py.yml"
    input:
        tsv="../resources/samples/samples.bed.grouped.count.treatment.tsv",
    output:
        out="../resources/samples/samples.bed.grouped.count.treatment.phrase.tsv",
    script:
        "../scripts/add_sample_type.py"


rule plot_read_counts:
    conda:
        "../envs/R.yml"
    input:
        tsv="../resources/samples/samples.bed.grouped.count.treatment.phrase.tsv",
    output:
        pdf='output/plots/read.counts.pdf'
    script:'../scripts/plots/plot.read.counts.R'


# rule merge_sample_treatments_all:
#     input:
#         expand(
#             'output/bedData/basepairConversionFreqSampleLabels/PBEH2-{sample_id}_Plasmid_{plasmid}_Strand_{strand}_call_{call_type}.conversion.freq.label.tsv',
#             zip, sample_id=BED_SAMPLES.sample_id, plasmid=BED_SAMPLES.plasmid, strand=BED_SAMPLES.strand, call_type=BED_SAMPLES.call_type
#         )
#     output:
#         'output/bedData/done.1.txt'
#     shell:'''
#     touch {output}
#     '''
