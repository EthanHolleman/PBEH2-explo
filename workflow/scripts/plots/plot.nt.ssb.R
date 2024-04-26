library(ggpubr)
library(RColorBrewer)
library(ggplot2)


df <- read.table(snakemake@input[['tsv']], sep='\t', header=TRUE)


# In this case want to break down plots by plasmid and compare the co
# vs post ssb samples and the controls on one plot. Seperate by strand
# and by call type

plot.freq <- function(df){

    title <- unique(df$plasmid)[1]
    p <- ggplot(df, aes(x=position, y=peak_freq, color=treatment_phrase)) +
            geom_line() +
            facet_grid(rows=vars(call_type), cols=vars(strand)) +
            theme_pubr() +
            scale_color_brewer(palette='Dark2') +
            xlim(0, 2000) +
            labs(title=title, x='position (bp)', y='R-loop freq')
    
    return(p)
}

plot.box <- function(df){

    # just compare the post vs co-SSB samples and 
    # limit scope to where most of signal actually is

    df <- subset(df, peak_freq > 0.001)

    df.ssb <- subset(
        df, treatment_phrase=='Nicked Linear Co-SSB Txn+' | treatment_phrase=='Nicked Linear Post-SSB Txn+'
        )
    
    df.ssb.pos <- subset(df.ssb, strand=='Pos')
    df.ssb.neg <- subset(df.ssb, strand=='Neg')

    plot <- function(data, title){

        ggplot(data, aes(x=plasmid, y=peak_freq, fill=treatment_phrase)) +
              geom_boxplot() + 
              theme_pubr() +
              scale_color_brewer(palette='Dark2') +
              labs(title=title, x='', y='R-loop frequency') +
              facet_wrap(~call_type) +
              coord_flip()
    }

    print(df.ssb.pos$call_type)

    p.pos <- plot(df.ssb.pos, 'Nt.BspQI SSB + strand')
    p.neg <- plot(df.ssb.neg, 'Nt.BspQI SSB - strand')

    print(p.pos)
    print(p.neg)

}


pdf(snakemake@output[['out']], width=12, height=6)

plasmids <- unique(df$plasmid)
print(plasmids)
for (each_plasmid in plasmids){

    df.plas <- subset(df, plasmid==each_plasmid)
    p <- plot.freq(df.plas)
    print(p)
}

plot.box(df)

dev.off()


