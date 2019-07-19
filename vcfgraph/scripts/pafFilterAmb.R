## Usage: Rscript pafFilterAmb.R al.paf graph.pdf chrName contigsToKeep.txt

args = commandArgs(TRUE)
## args = c('minimap2/HG02723_hg38.paf',
##          'test.pdf', 'chr10', 'test.txt')

library(dplyr)
library(ggplot2)
library(GenomicRanges)

## Read PAF file
paf = read.table(args[1], as.is=TRUE, sep='\t', fill=TRUE,
                 row.names=NULL,
                 colClasses=c('character', rep('integer', 3),
                              rep('character', 2), rep('integer', 6),
                              rep('NULL', 13)))
colnames(paf) = c('qseqn','qlen', 'qstart','qend', 'qstand', 'tseqn','tlen', 'tstart','tend', 'nres', 'allen', 'mapq')

## Remove alternate ref contigs
paf = subset(paf, !grepl('_', tseqn))

## Remove alignment with mapq=0 and make genomic ranges
pafq.gr = with(subset(paf, mapq>0), GRanges(qseqn, IRanges(qstart, qend), qlen=qlen, tseqn=tseqn))

## Split by tseqn and reduce
pafq.l = split(pafq.gr, pafq.gr$tseqn)
pafq.l = reduce(pafq.l)
pafq.l = lapply(names(pafq.l), function(tseqn){
  gr = pafq.l[[tseqn]]
  gr$tseqn = tseqn
  gr
})
pafq.r = do.call(c, pafq.l)

## Add query length back
ql = paf %>% select(qseqn, qlen) %>% unique
qlConv = ql$qlen
names(qlConv) = ql$qseqn
pafq.r$qlen = qlConv[as.character(seqnames(pafq.r))]

## Compute coverage metrics
paf.q = pafq.r %>% as.data.frame %>%
  group_by(seqnames, tseqn, qlen) %>%
  summarize(width=sum(width), prop=width/qlen[1]) %>%
  group_by(seqnames) %>% mutate(rank=rank(prop))

paf.top2 = paf.q %>% group_by(seqnames, qlen) %>%
  arrange(desc(prop)) %>% do(head(.,2)) %>%
  summarize(prop1=prop[1], prop2=prop[2],
            tseqn1=tseqn[1], tseqn2=tseqn[2]) %>%
  mutate(prop2=ifelse(is.na(prop2), 0, prop2))

## Graphs
pdf(args[2], 9, 7)
ggplot(paf, aes(mapq)) + geom_histogram() + theme_bw() +
  ylab('alignment')
paf %>% group_by(qseqn) %>% summarize(mapq=max(mapq)) %>%
  ggplot(aes(mapq)) + geom_histogram() + theme_bw() +
  xlab('maximum mapq') + ylab('contig')
paf.top2 %>% mutate(qlen.class=cut(qlen, c(1,1000,1e4,1e5,1e6,Inf))) %>%
  ggplot(aes(x=prop1, y=prop2, color=qlen.class)) +
  geom_point(alpha=.7) +
  scale_color_brewer(name='contig length', palette='Set1') + 
  theme_bw() + xlab('proportion aligning to best ref chr') +
  ylab('proportion aligning to second best ref chr') +
  scale_x_continuous(breaks=seq(0,1,.2)) +
  scale_y_continuous(breaks=seq(0,1,.2))
ggplot(paf.top2, aes(x=qlen, y=prop1-prop2, colour=prop1)) +
  geom_point(alpha=.7) + theme_bw() + scale_x_log10() +
  xlab('contig length') + ylab('delta proportion aligning to top 2 ref chr') +
  scale_color_gradient(low='black', name='proportion\naligning\nto best ref chr')
## Chromosome of interest
paf.top2 %>% mutate(qlen.class=cut(qlen, c(1,1000,1e4,1e5,1e6,Inf))) %>%
  filter(tseqn1==args[3]) %>% 
  ggplot(aes(x=prop1, y=prop2, color=qlen.class)) +
  geom_point(alpha=.7) +
  scale_color_brewer(name='contig length', palette='Set1') + 
  theme_bw() + xlab('proportion aligning to best ref chr') +
  ylab('proportion aligning to second best ref chr') +
  scale_x_continuous(breaks=seq(0,1,.2)) +
  scale_y_continuous(breaks=seq(0,1,.2)) +
  ggtitle(args[3])
paf.top2 %>% filter(tseqn1==args[3]) %>% 
  ggplot(aes(x=qlen, y=prop1-prop2, colour=prop1)) +
  geom_point(alpha=.7) + theme_bw() + scale_x_log10() +
  xlab('contig length') + ylab('delta proportion aligning to top 2 ref chr') +
  scale_color_gradient(low='black', name='proportion\naligning\nto best ref chr') + 
  ggtitle(args[3])
dev.off()

## Contigs to keep
cont.sel = paf.top2 %>% filter(tseqn1==args[3]) %>% .$seqnames %>% sort
write(as.numeric(levels(cont.sel)[cont.sel]), file=args[4], ncolumns=1)
