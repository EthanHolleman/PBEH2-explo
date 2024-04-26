library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(tidyr)


# plotting functions

plot.read.counts <- function(df.exp, plot_title, y_val){

    p <- ggplot(df.exp, aes_string(x='plasmid', fill='strand', y=y_val)) +
        geom_bar(stat='identity', position='dodge') +
        theme_pubr() +
        labs(title=plot_title, x='Plasmid', y='Total reads') +
        facet_wrap(~treatment_phrase) +
        coord_flip() +
        scale_fill_brewer(palette='Dark2')
    
    p

}


# do the plotting here 

df <- read.table(snakemake@input[['tsv']], sep='\t', header=TRUE)
experiments <- unique(df$Experiment)

pdf(snakemake@output[['pdf']], width=10, height=6)

for (each_exp in experiments){

    df.exp <- subset(df, Experiment==each_exp)
    df.exp <- subset(df.exp, call_type=='high')  # read counts for call types are same


    df.exp$by.plasmid <- df.exp$read_count / df.exp$Total.plasmids.in.sample
    title=unique(df.exp$Experiment_name)[1]


    p1 <- plot.read.counts(df.exp, title, 'read_count')
    #p2 <- plot.read.counts(df.exp, title, 'by.plasmid')

    #p2 <- labs(subtitle='Per plasmid')

    print(p1)
    #print(p2)

}

dev.off()













