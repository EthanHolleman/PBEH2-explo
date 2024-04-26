library(ggpubr)
library(RColorBrewer)
library(ggplot2)


df <- read.table(snakemake@input[['tsv']], sep='\t', header=TRUE)

# should only be 1 plasmid in these samples so don't need to break down
# plots by plasmid. Here we most want to compare between linearization
# before and after transcription between the supercoiled and linear sample
# types

df.pos <- subset(df, strand=='Pos')
df.neg <- subset(df, strand=='Neg')

df.pos.linear <- subset(df.pos, Linear!='-')
df.pos.sc <- subset(df.pos, Linear=='-')

df.neg.linear <- subset(df.neg, Linear!='-')
df.neg.sc <- subset(df.neg, Linear=='-')


freq.plot <- function(df.strand, title, sub.title){


    ggplot(df.strand, aes(x=position, y=peak_freq, color=treatment_phrase)) +
           geom_line() +
           theme_pubr() +
           scale_color_brewer(palette='Dark2') +
           xlim(0, 1400) +
           labs(x='position (bp)', y='R-loop frequency', title=title, subtitle=sub.title) +
           facet_wrap(~call_type, ncol=1) + ylim(0, 1)

}

pdf(snakemake@output[['out']], width=13, height=6)

print(freq.plot(df.pos.linear, 'SSB stabilization + strand', 'Linear'))
print(freq.plot(df.pos.sc, 'SSB stabilization + strand', 'SC'))
print(freq.plot(df.neg.linear, 'SSB stabilization - strand', 'Linear'))
print(freq.plot(df.neg.sc, 'SSB stabilization - strand', 'SC'))

dev.off()
