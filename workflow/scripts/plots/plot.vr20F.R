library(ggplot2)
library(ggpubr)
library(RColorBrewer)

df <- read.table(snakemake@input[['tsv']], sep='\t', header=TRUE)

VR.START <- 349
VR.END <- VR.START + 200

df.plus <- subset(df, strand=='Pos')
df.neg <- subset(df, strand=='Neg')

df.plus.vr <- subset(df.plus, position >= VR.START & position <= VR.END)
df.neg.vr <- subset(df.plus, position >= VR.START & position <= VR.END)

df.plus.vr.SC <- subset(
        df.plus.vr, treatment_phrase=='SC Txn+' | treatment_phrase=='SC Txn-' | treatment_phrase=='SC RnaseH'
        )

print(nrow(df.plus.vr.SC))
print(nrow(df.plus.vr))
print(nrow(df.plus))


pdf(snakemake@output[['out']], width = 12, height = 16)

# Frequency type position plots
# =============================================================================

p.plus <- ggplot(df.plus, aes(x=position, y=peak_freq, color=plasmid)) +
        facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') +
        labs(title='Conversion frequency + strand', x='', y='R-loop freq', color='') +
        geom_line() +
        ylim(0, 0.2) +
        xlim(0, 2000) +
        geom_vline(xintercept = VR.START,  size=1, linetype = "dashed") +
        geom_vline(xintercept = VR.END, size=1, linetype = "dashed") +
        labs(title='R-loop frequency VR20F samples + strand')

p.neg <- ggplot(df.neg, aes(x=position, y=peak_freq, color=plasmid)) +
        geom_line() +
        facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
        theme_pubr() +
        scale_color_brewer(palette='Dark2') + 
        labs(title='Conversion frequency - strand', x='', y='R-loop freq', color='') +
        ylim(0, 0.2) +
        xlim(0, 2000) + 
        geom_vline(xintercept = VR.START, size=0.5, linetype = "dashed") +
        geom_vline(xintercept = VR.END, size=0.5, linetype = "dashed") +
        labs(title='R-loop frequency VR20F samples - strand')

p.plus.zoom <- p.plus + xlim(VR.START+1, VR.END-1) + 
                labs(title='R-loop frequency VR20F VR region + strand')
p.neg.zoom <- p.neg + xlim(VR.START+1, VR.END-1) + 
                labs(title='R-loop frequency VR20F VR region - strand')

# Boxplots comparing R-loop frequency within the VR region specifically
# =============================================================================

p.plus.box <- ggplot(df.plus.vr, aes(x=class, fill=as.factor(version), y=peak_freq)) +
         geom_boxplot() +
         scale_fill_brewer(palette='Dark2') +
         theme_pubr() +
         facet_grid(rows=vars(treatment_phrase), cols=vars(call_type))


p.neg.box <- ggplot(df.neg.vr, aes(x=class, fill=as.factor(version), y=peak_freq)) +
         geom_boxplot() +
         scale_fill_brewer(palette='Dark2') +
         theme_pubr() +
         facet_grid(rows=vars(treatment_phrase), cols=vars(call_type))

print(p.plus)
print(p.neg)
print(p.plus.zoom)
print(p.neg.zoom)
print(p.plus.box)
print(p.neg.box)

dev.off()

# Focusing in on specific samples 
# =============================================================================

comparisons <- list(
        c(1, 2),
        c(2, 3),
        c(1, 3)
)

p.sc.box <- ggplot(df.plus.vr.SC, aes(x=class, fill=as.factor(version), y=peak_freq)) +
            geom_boxplot() + 
            stat_summary(fun.y=median, geom="point", size=2, color="black") +
            scale_fill_brewer(palette='Dark2') +
            stat_compare_means(comparisons = comparisons, method='t.test') +
            theme_pubr() +
            facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
            labs(y='VR R-loop frequency SC samples VR20F SC samples', fill='Edition', x='Folding class')

# select only positions at which at least one peak was measured (peak length)
# must be greater than 0

df.plus.vr.SC.len <- subset(df.plus.vr.SC, mean_len > 1)

p.sc.box.len <- ggplot(df.plus.vr.SC.len, aes(x=class, fill=as.factor(version), y=mean_len)) +
            geom_boxplot() + 
            stat_summary(fun.y=median, geom="point", size=2, color="black") +
            scale_fill_brewer(palette='Dark2') +
            stat_compare_means(comparisons = comparisons, method='t.test') +
            theme_pubr() +
            facet_grid(rows=vars(treatment_phrase), cols=vars(call_type)) +
            labs(y='VR R-loop R-loop length VR20F SC samples', fill='Edition', x='Folding class')


ggsave('output/plots/VR20F.SC.TXN+.Box.png', p.sc.box, dpi=400, width=10, height=10)
ggsave('output/plots/VR20F.SC.TXN+.Box.LENGTH.png', p.sc.box.len, dpi=400, width=10, height=10)


