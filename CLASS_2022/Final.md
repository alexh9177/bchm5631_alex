Genome-wide DNA binding in HEPG2
================
Alex Hirano
4/17/2022

# Goal:

Here we aim to download all available DNA binding protein (DBP) profiles
in a single cell state (measured by ChIP-seq) . This will allow us to
investigate the binding properties of hundreds of DBPs in the same
cellular context or background. We aim to address several questions: (i)
What are the number of peaks and genome coverage for each DBP? (ii) What
are the binding preferences for promoters, gene-bodies and intergenic
genomic regions? (iii) What are the similarities and differences across
DBPs based on their genome-wide binding profiles genome-wide? (iv) What
properties or preferences do promoters have for binding events. (iv) Are
there reservoir promoters in HepG2 as defined in k562 previously? (v)
How does binding to a promoter affect the transcriptional output of that
promoter?

To address these questions we have curated a set of X,000 ChIPs-eq data
sets comprised of 486 DBPs in HEPG2 cells from the ENCODE consortrium.
We required duplicate ChIP-seq experiments for a given DBP and other
criterion that can be found here :

<https://www.encodeproject.org/report/?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>

## These samples were selected on the following criteria:

1.  “chromatin” interaction data, then DNA binding data, cell line
    HEPG2, “TF-Chip-seq”.
2.  We further selected “TF Chip-seq”, “Control chip-seq” and “Histone
    Chip-seq”.
3.  We selected several read lengths to get the most DNA binding
    proteins (DBPs)
4.  Read lengths: 100, 76, 75, 36
5.  ONLY SINGLE END READS (this eliminates 54 samples)

### Experimental data was downloading by (ENCODE report.tsv):

<https://www.encodeproject.org/report.tsv?type=Experiment&status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended>

### The FASTQ files were downloaded with:

“<https://www.encodeproject.org/metadata/?status=released&assay_slims=DNA+binding&biosample_ontology.term_name=HepG2&assay_title=TF+ChIP-seq&biosample_ontology.classification=cell+line&files.read_length=100&files.read_length=76&files.read_length=75&files.read_length=36&assay_title=Control+ChIP-seq&assay_title=Histone+ChIP-seq&files.run_type=single-ended&type=Experiment>”

MD5sums were checked with all passing (see encode\_file\_info function
to reterive MD5Sum values that are not available from the encode portal
(/util)

### Processing data:

We processed all the read alignments and peak calling using the NF\_CORE
ChIP-seq pipeline: (nfcore/chipseq v1.2.1)

## Next we created consensus peaks that overlap in both replicates

Our strategy was to take peaks in each replicate and find all
overlapping peak windows. We then took the union length of the
overlapping range in each peak window.

``` r
# create_consensus_peaks requires an annotation .GTF file - loading in Gencode v32 annotations.
gencode_gr <- rtracklayer::import("/scratch/Shares/rinnclass/CLASS_2022/data/genomes/gencode.v32.annotation.gtf")
# Creating consensus peaks function to create a .bed file of overlapping peaks in each replicate.
# /util/intersect_functions.R
# TODO run this only on final knit
# create_consensus_peaks <- create_consensus_peaks(broadpeakfilepath = "/scratch/Shares/rinnclass/CLASS_2022/data/test_work/all_peak_files")
# exporting consensus peaks .bed files
# TODO run this only on final knit
# for(i in 1:length(consensus_peaks)) {
#  rtracklayer::export(consensus_peaks[[i]],
#                     paste0("/scratch/Shares/rinnclass/CLASS_2022/JR/CLASS_2022/class_exeRcises/analysis/11_consensus_peaks/consensus_peaks/",
#                             names(consensus_peaks)[i],
#                             "_consensus_peaks.bed"))
# }
```

# loading in consensus peaks to prevent rerunning create\_consensus\_peaks function

``` r
# Loading in files via listing and rtracklayer import
consensus_fl <- list.files("/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/class_exeRcises/analysis/11_consensus_peaks/consensus_peaks", full.names = T)
# importing (takes ~5min)
consensus_peaks <- lapply(consensus_fl, rtracklayer::import)
# cleaning up file names
names(consensus_peaks) <- gsub("/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/class_exeRcises/analysis/11_consensus_peaks/consensus_peaks/|_consensus_peaks.bed","", consensus_fl)
# Filtering consensus peaks to those DBPs with at least 250 peaks
num_peaks_threshold <- 250
num_peaks <- sapply(consensus_peaks, length)
filtered_consensus_peaks <- consensus_peaks[num_peaks > num_peaks_threshold]
# Result: these were the DBPs that were filtered out.
filtered_dbps <- consensus_peaks[num_peaks < num_peaks_threshold]
names(filtered_dbps)
```

    ##  [1] "CEBPZ"    "GPBP1L1"  "H3K27me3" "HMGA1"    "IRF3"     "MLLT10"  
    ##  [7] "MYBL2"    "NCOA5"    "RNF219"   "RORA"     "ZBTB3"    "ZFP36"   
    ## [13] "ZFP62"    "ZMAT5"    "ZNF10"    "ZNF17"    "ZNF260"   "ZNF382"  
    ## [19] "ZNF48"    "ZNF484"   "ZNF577"   "ZNF597"   "ZNF7"

``` r
# We have this many remaining DBPs
length(filtered_consensus_peaks)
```

    ## [1] 460

## Now we will determine the peak number and genome coverage for each DBP.

``` r
# Let's start with loading in the number of peaks each DBP has -- using length.
num_peaks_df <- data.frame("dbp" = names(filtered_consensus_peaks),
                           "num_peaks" = sapply(filtered_consensus_peaks, length))
# total genomic coverage of peaks for each dbp
num_peaks_df$total_peak_length <- sapply(filtered_consensus_peaks, function(x) sum(width(x)))
# Plotting distribution of peak number per dbp
hist(num_peaks_df$total_peak_length)
```

![](Final_files/figure-gfm/peak%20number%20and%20coverage%20per%20DBP-1.png)<!-- -->

# Now we will create promoter annotations for lncRNA and mRNA and both.

We have created a function get\_promter\_regions that has up and
downstream parameters

``` r
# creating lncRNA and mRNA promoters
lncrna_mrna_promoters <- get_promoter_regions(gencode_gr, biotype = c("lncRNA", "protein_coding", upstream = 3000, downstream = 3000))
names(lncrna_mrna_promoters) <- lncrna_mrna_promoters$gene_id
rtracklayer::export(lncrna_mrna_promoters, "analysis/results/lncRNA_mrna_promoters.gtf")
# creating lncRNAs promoter
lncrna_promoters <- get_promoter_regions(gencode_gr, biotype = "lncRNA", upstream = 3000, downstream = 3000) 
names(lncrna_promoters) <- lncrna_promoters$gene_id
rtracklayer::export(lncrna_promoters, "analysis/results/lncRNA_promoters.gtf")
# creating mRNA promoters
mrna_promoters <- get_promoter_regions(gencode_gr, biotype = "protein_coding", upstream = 3000, downstream = 3000)
names(mrna_promoters) <- mrna_promoters$gene_id
rtracklayer::export(lncrna_promoters, "analysis/results/mRNA_promoters.gtf")
# creating all genebody annotation
lncrna_mrna_genebody <- gencode_gr[gencode_gr$type == "gene" & 
                                     gencode_gr$gene_type %in% c("lncRNA", "protein_coding")]
names(lncrna_mrna_genebody) <- lncrna_mrna_genebody$gene_id
rtracklayer::export(lncrna_mrna_genebody, "analysis/results/lncrna_mrna_genebody.gtf")
# creating lncRNA genebody annotation
lncrna_genebody <- gencode_gr[gencode_gr$type == "gene" & 
                                gencode_gr$gene_type %in% c("lncRNA")]
names(lncrna_genebody) <- lncrna_genebody$gene_id
rtracklayer::export(lncrna_mrna_genebody, "analysis/results/lncrna_genebody.gtf")
# creating mRNA genebody annotation
mrna_genebody <- gencode_gr[gencode_gr$type == "gene" & 
                              gencode_gr$gene_type %in% c("protein_coding")]
names(mrna_genebody) <-mrna_genebody$gene_id
rtracklayer::export(lncrna_mrna_genebody, "analysis/results/mrna_genebody.gtf")
```

# Determining the overlaps of chip peaks with promoters and genebodys

``` r
# creating index to subset lncRNA and mRNA annotations
lncrna_gene_ids <- lncrna_mrna_genebody$gene_id[lncrna_mrna_genebody$gene_type == "lncRNA"]
mrna_gene_ids <- lncrna_mrna_genebody$gene_id[lncrna_mrna_genebody$gene_type == "protein_coding"]
# using count peaks per feature returns number of annotation overlaps for a given DBP (takes ~5min)
promoter_peak_counts <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_peaks, type = "counts")
# adding data to num_peaks_df
num_peaks_df$peaks_overlapping_promoters <- rowSums(promoter_peak_counts)
num_peaks_df$peaks_overlapping_lncrna_promoters <- rowSums(promoter_peak_counts[,lncrna_gene_ids])
num_peaks_df$peaks_overlapping_mrna_promoters <- rowSums(promoter_peak_counts[,mrna_gene_ids])
# gene body overlaps 
genebody_peak_counts <- count_peaks_per_feature(lncrna_mrna_genebody, 
                                                filtered_consensus_peaks, 
                                                type = "counts")
# adding data to num_peaks_df
num_peaks_df$peaks_overlapping_genebody <- rowSums(genebody_peak_counts)
num_peaks_df$peaks_overlapping_lncrna_genebody <- rowSums(genebody_peak_counts[,lncrna_gene_ids])
num_peaks_df$peaks_overlapping_mrna_genebody <- rowSums(genebody_peak_counts[,mrna_gene_ids])
write_csv(num_peaks_df, "analysis/results/num_peaks_df.csv")
```

# Plotting peak annotation features for DBPs

``` r
num_peaks_df <- read_csv("analysis/results/num_peaks_df.csv")
# Distribution of peak numbers of all 460 DBPs
ggplot(num_peaks_df, aes(x = num_peaks)) + 
  geom_histogram(bins = 70)
```

![](Final_files/figure-gfm/plotting%20peak%20annotation%20features-1.png)<!-- -->

``` r
# Plotting number of peaks versus total genome coverage
ggplot(num_peaks_df, aes(x = num_peaks, y = total_peak_length)) +
  geom_point() + 
  geom_smooth(method = "gam", se = TRUE, color = "black", lty = 2)+
         
  ylab("BP covered") +
  xlab("Number of peaks") +
  ggtitle("Peak count vs. total bases covered")
```

![](Final_files/figure-gfm/plotting%20peak%20annotation%20features-2.png)<!-- -->

``` r
ggsave("analysis/figures/peak_num_vs_coverage.pdf")
# Plotting number of peaks versus peaks overlapping promoters
ggplot(num_peaks_df,
       aes(x = num_peaks, y = peaks_overlapping_promoters)) +
  xlab("Peaks per DBP") +
  ylab("Number of peaks overlapping promoters") +
  ggtitle("Relationship Between Number of DBP Peaks and Promoter Overlaps")+
  geom_point() +
  geom_abline(slope = 1, linetype="dashed") +
  geom_smooth(method = "lm", se=FALSE, formula = 'y ~ x',
              color = "#a8404c") +
  stat_regline_equation(label.x = 35000, label.y = 18000) +
  ylim(0,60100) +
  xlim(0,60100)
```

![](Final_files/figure-gfm/plotting%20peak%20annotation%20features-3.png)<!-- -->

``` r
ggsave("analysis/figures/3_peak_num_vs_promoter_coverage.pdf")
# Plotting peak overlaps with genebody
ggplot(num_peaks_df,
       aes(x = num_peaks, y = peaks_overlapping_genebody)) +
  xlab("Peaks per DBP") +
  ylab("Number of peaks overlapping genes") +
  ggtitle("Relationship Between Number of DBP Peaks and Gene Body Overlaps")+
  geom_point() +
  geom_abline(slope = 1, linetype="dashed") +
  geom_smooth(method = "lm", se=F, formula = 'y ~ x',
              color = "#a8404c") +
  stat_regline_equation(label.x = 35000, label.y = 18000) +
  ylim(0,60100) +
  xlim(0,60100)
```

![](Final_files/figure-gfm/plotting%20peak%20annotation%20features-4.png)<!-- -->

``` r
ggsave("analysis/figures/4_peak_num_vs_gene_body_coverage.pdf")
# there is a large amount of data explained (almost all by genebodys)
# Let's see what percentage of the genome genebodys cover:
reduced_gene_bodies <- gencode_gr[gencode_gr$type == "gene"] %>%
  GenomicRanges::reduce() %>%
  width() %>%
  sum()
# percentage of gene bodies in genome
reduced_gene_bodies/3.2e9
```

    ## [1] 0.589159

# Counting the number of overlaps at each promoter

Promoters are the cols and DBPs rows thus we can retrieve the number of
binding events at each promoter unlike the “counts parameter” that just
gives total number of overlaps

``` r
# Creating matrix of promoters(annotation feature) as cols and DBPs as rows (takes ~5min)
promoter_peak_occurence <- count_peaks_per_feature(lncrna_mrna_promoters, filtered_consensus_peaks, 
                                               type = "occurrence")
# test to make sure everything is in right order
stopifnot(all(colnames(promoter_peak_occurence) == lncrna_mrna_promoters$gene_id))
# Formatting final data.frame from peak occurrence matrix
peak_occurence_df <- data.frame("gene_id" = colnames(promoter_peak_occurence),
                                "gene_name" = lncrna_mrna_promoters$gene_name,
                                "gene_type" = lncrna_mrna_promoters$gene_type,
                                "chr" = lncrna_mrna_promoters@seqnames,   
                                "1.5_kb_up_tss_start" = lncrna_mrna_promoters@ranges@start,
                                "strand" = lncrna_mrna_promoters@strand,
                                "number_of_dbp" = colSums(promoter_peak_occurence))
# exporting
#write_csv(peak_occurence_df, "analysis/results/peak_occurence_dataframe.csv")
```

# WHERE DO YOU WANT TO GO FROM HERE ?

# lncRNA and mRNA clustering comparison

We wanted to examine the difference between lncRNA and mRNA clustering
w/respect to POLR2A and POLR2A S5/S2

POLR2A codes for RNA polymerase II subunit 1, where S2 and S5 indicate
phosphorylations on serines 2 and 5 of the hepapeptide repeats in the
C-terminal domain CTD.

``` r
# Used to build the peak occurrence matrix
write.table(promoter_peak_occurence, "/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/analysis/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")

# Build peak occurrence matrix
promoter_peak_occurence_matrix = read.table("analysis/results/lncrna_mrna_promoter_peak_occurence_matrix.tsv")
promoter_peak_occurence_matrix = as.matrix(promoter_peak_occurence_matrix)
peak_occurence_matrix <- promoter_peak_occurence_matrix[rowSums(promoter_peak_occurence_matrix) > 250, ]

# Import the peak occurrence df then combine with the peak occurrence matrix for lncRNA
peak_occurence_df <- read.csv("/scratch/Shares/rinnclass/CLASS_2022/alhi9177/bchm5631_alex/CLASS_2022/analysis/results/peak_occurence_dataframe.csv")
lncrna_promoters <- rtracklayer::import("analysis/results/lncRNA_promoters.gtf")
mrna_promoters <- rtracklayer::import("analysis/results/mRNA_promoters.gtf")
lncrna_peak_occurence <- peak_occurence_matrix[,lncrna_promoters$gene_id]
bin_hier_lncrna <- hclust(dist(lncrna_peak_occurence, method = "binary"))
ggdendro::ggdendrogram(bin_hier_lncrna, rotate = TRUE, size = 10)
```

![](Final_files/figure-gfm/lncRNA%20promoter%20clustering-1.png)<!-- -->

``` r
#ggsave("analysis/figures/lncrna_hclust_binary_dist.pdf", height = 49, width = 6)

# Repeat the same for mRNA
mrna_peak_occurence <- peak_occurence_matrix[,mrna_promoters$gene_id]
bin_hier_mrna <- hclust(dist(mrna_peak_occurence, method = "binary"))
ggdendro::ggdendrogram(bin_hier_mrna, rotate = TRUE,  size = 10)
```

![](Final_files/figure-gfm/lncRNA%20promoter%20clustering-2.png)<!-- -->

``` r
#ggsave("analysis/figures/mrna_hclust_binary_dist.pdf", height = 44, width = 6)

knitr::include_graphics('analysis/figures/lncrna_hclust_binary_dist.pdf')
```

![](analysis/figures/lncrna_hclust_binary_dist.pdf)<!-- -->

``` r
knitr::include_graphics('analysis/figures/mrna_hclust_binary_dist.pdf')
```

![](analysis/figures/mrna_hclust_binary_dist.pdf)<!-- -->

# Results

There doesn’t seem to be a difference in preference between POLR2A
phospho S2/S5 binding between mRNA and lncRNA in Hep32. For both plots,
POLR2A S2 is flanked by ASR2L, FUBP1, and TRAFD1 above; and POL2A, TP53,
and FOXJ3 below. POLR2A S5 is flanked by ZFN652, ZFN629, and ZFN 710
above and PAF1, ZBTB7A, and NFE2L1 below. Both of the POL2RA genes
diverge at the same location between mRNA and lncRNA, suggesting no
preference in binding affinity for either.

# Dendrogram plot

We wanted to explore sample similarity after generating bin\_hier
(cross-correlation across all samples)

``` r
# Generate peak occurrence and bin_hier
peak_occurence_matrix[1:2,1:5]
```

    ##      ENSG00000243485.5 ENSG00000237613.2 ENSG00000186092.6 ENSG00000238009.6
    ## ADNP                 0                 0                 0                 0
    ## AFF4                 0                 0                 0                 0
    ##      ENSG00000239945.1
    ## ADNP                 0
    ## AFF4                 0

``` r
peak_occurence_dist <- dist(peak_occurence_matrix, method = "binary")
peak_occurence_dist[1:10]
```

    ##  [1] 0.8088089 0.8358010 0.8238244 0.8349438 0.9438168 0.8293561 0.8316339
    ##  [8] 0.8396660 0.8504692 0.8194388

``` r
# Generating a ggdendrogram
bin_hier <- hclust(peak_occurence_dist, method = "complete")

ggdendro::ggdendrogram(bin_hier, rotate = FALSE,  size = 3, 
                       theme_dendro = TRUE) +
   coord_flip() +
   scale_y_continuous() +
   scale_x_continuous(position = "top") +
   scale_x_continuous(breaks = seq_along(bin_hier$labels[bin_hier$order]),
             labels = bin_hier$labels[bin_hier$order], position = "top",
             expand = c(0,0)) +
   theme(axis.text.x = element_text(angle = 90, hjust  = 1)) + 
   theme(axis.text.y = element_text(angle = 0,hjust = 1)) +
   scale_y_reverse(expand = c(0.01, 0)) +
   theme(
     plot.background = element_blank(),
     panel.grid.major = element_blank(),
   panel.grid.minor = element_blank(),
     panel.border = element_blank()
   )
```

![](Final_files/figure-gfm/Generating%20dendrogram%20plot-1.png)<!-- -->

``` r
#ggsave("analysis/figures/ggdendro_plot.pdf", height = 50, width = 12, limitsize = F)
knitr::include_graphics('analysis/figures/ggdendro_plot.pdf')
```

![](analysis/figures/ggdendro_plot.pdf)<!-- -->

# Results

This generates a complete dendrogram but there’s not much we can
elucidate from this figure as-is. Subsequent figures will cover a
clustering comparison between lncRNA and mRNA.

# Making a heat map of high-binding promoters

We wanted to examine the occurrence of high-binding promoters for our
3kb window using a heat map

``` r
# Generating the high binders
high_binders <- promoter_peak_occurence_matrix[,colSums(promoter_peak_occurence_matrix) > 300]
peak_occurence_matrix_test <- high_binders[rowSums(promoter_peak_occurence_matrix) > 250, ]
ncol(high_binders)
```

    ## [1] 5663

``` r
min(rowSums(peak_occurence_matrix_test)) 
```

    ## [1] 26

``` r
which.min(rowSums(high_binders))
```

    ## H3K9me3 
    ##      94

``` r
# Which DBP binds the high binders? The low binders?
max(rowSums(high_binders))
```

    ## [1] 5663

``` r
table(max(rowSums(high_binders)))
```

    ## 
    ## 5663 
    ##    1

``` r
which.max(rowSums(high_binders))
```

    ## ARID4B 
    ##      9

``` r
# Generate and save heatmap
xx <- pheatmap(high_binders, show_colnames = FALSE, clustering_distance_rows = "binary", clustering_distance_cols = "binary")
```

![](Final_files/figure-gfm/Heat%20map%20clustering%20high-binding%20promoters-1.png)<!-- -->

``` r
save_pheatmap_pdf <- function(x, filename, width=20, height=20) {
   stopifnot(!missing(x))
   stopifnot(!missing(filename))
   pdf(filename, width=width, height=height)
   grid::grid.newpage()
   grid::grid.draw(x$gtable)
   dev.off()
}
#save_pheatmap_pdf(xx, "analysis/figures/test.pdf")
knitr::include_graphics('analysis/figures/test.pdf')
```

![](analysis/figures/test.pdf)<!-- -->

# Results

These results make sense - filtering for high binders and binary, you’re
either highly bound (1, red) or not highly bound (0, blue). The majority
of the promoters here seem to be highly bound, especially towards the
middle of the dendrogram. Numerous promoters seem to be high-binding for
HepG2.

POLR2AphosS5 seems to be more high-binding than POLR2AphosS2, which
makes sense as the phosS2 segment (largest segment of RNA Pol II) is
actively involved in transcription and not constantly bound.
