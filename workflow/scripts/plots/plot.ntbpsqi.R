library(ggplot2)
library(ggpubr)
library(RColorBrewer)


df <- read.table(snakemake@input[['tsv']], sep='\t', header=TRUE)

df.plus <- subset(df, strand=='Pos')
df.neg <- subset(df, strand=='Neg')


pdf(snakemake@output[['out']], width = 12, height = 8)

# Peak frequency plots

p.plus <- ggplot(df.plus, aes(x=position, y=peak_freq, color=treatment_phrase)) +
        facet_grid(rows=vars(plasmid), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') +
        labs(title='Conversion frequency + strand', x='', y='R-loop freq', color='') +
        geom_line() +
        ylim(0, 0.2) +
        xlim(0, 2000)

p.neg <- ggplot(df.neg, aes(x=position, y=peak_freq, color=treatment_phrase)) +
        geom_line() +
        facet_grid(rows=vars(plasmid), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') + 
        labs(title='Conversion frequency - strand', x='', y='R-loop freq', color='') +
        ylim(0, 0.2) +
        xlim(0, 2000)

# Peak length plots

p.plus.len <- ggplot(df.plus, aes(x=position, y=mean_len, color=plasmid)) +
        facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') +
        labs(title='Mean peak length + strand', x='', y='Mean peak length', color='') +
        geom_line() +
        xlim(0, 2000)

p.neg.len <- ggplot(df.neg, aes(x=position, y=mean_len, color=plasmid)) +
        geom_line() +
        facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') + 
        labs(title='Mean peak length - strand', x='', y='Mean peak length', color='') +
        xlim(0, 2000)

print(p.plus)
print(p.neg)

print(p.plus.len)
print(p.neg.len)

dev.off()





