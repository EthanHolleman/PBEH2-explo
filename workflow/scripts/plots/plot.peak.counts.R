library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(dplyr)


plot.peak.counts <- function(experiment.num, df){
    
    
    # count number of peaks by sample and divide by total plasmids in the sample
    df <- subset(df, Experiment==experiment.num)
    print(head(df))
    experiment.name <- unique(df$Experiment_name)
    total.plasmid <- unique(df$Total.plasmids.in.sample)
    count.df <- df %>%
              group_by(strand, call_type, sample_id, Txn) %>%
              summarise(count = n())
    count.df$plasmid.count <- count.df$count
    
    # make a tile plot of what the treatments of each of the samples are
    
    
    ggplot(count.df, aes(x=as.factor(sample_id), y=plasmid.count, fill=call_type)) + 
            geom_bar(position='dodge', stat='identity') + 
            scale_fill_brewer(palette='Dark2') +
            theme_pubr() + 
            labs(
                fill='Call type', y='Peak count per plasmid', 
                x='Sample ID', title=experiment.name
            ) +
            facet_wrap(~Txn)
}


peaks <- read.table(snakemake@input[['tsv']], header = TRUE, sep = "\t")
experiments <- unique(peaks$Experiment)

pdf(snakemake@output[['out']])

for (each_exp in experiments){

    df <- subset(peaks, Experiment==each_exp)
    p <- plot.peak.counts(each_exp, df)
    print(p)
}

dev.off()


