library("sleuth")
library("gridExtra")
library("cowplot")
library("biomaRt")
library("ggplot2")
library("dplyr")
library("ggrepel")

################################################################################
# Prepare the files required for the analysis & get the description of the transcripts from ensembl 

bench_isols <- c('21-0001', '21-0003', '21-0005', '21-0007', '21-0008', '21-0009', '21-0013')
amb_isols <- c('21-0006','21-0010','21-0016','21-0021')
kalm_isols <- c('21-0004','21-0014','21-0015','21-0018','21-0019','21-0020')

kal_dirs <- file.path("abundance", c(bench_isols,amb_isols))

s2c <- read.table(file.path("kallisto_samples.tsv"),
                  header = TRUE,
                  stringsAsFactors = FALSE,
                  sep = "\t") |>
  dplyr::filter(variety %in% c('benchmark','amboise'))

s2c <- dplyr::mutate(s2c, path = kal_dirs)

################################################################################
# list datasets from Ensembl and select the one of interest

marts <- listMarts(host = "https://fungi.ensembl.org")
marts
fungi_mart <- useMart("fungi_mart", host = "https://fungi.ensembl.org" )
listDatasets(fungi_mart)
fungi_mart <- useMart("fungi_mart", dataset = "pstriiformis_eg_gene", host="https://fungi.ensembl.org" )

# Now, view the attributes, and select the ones of interest.
# Definitely need the transcript_id because that's the one we used in Kallisto
listAttributes(fungi_mart)
t2g <- getBM(attributes = c("ensembl_transcript_id",
                            "ensembl_gene_id",
                            "chromosome_name",
                            "external_gene_name"),
             mart =fungi_mart)

ttg <- dplyr::rename(t2g, target_id= ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, chr_name = chromosome_name)

s2c$var <- as.factor(s2c$var)
s2c$variety <- factor(s2c$variety, levels = c('benchmark', 'amboise'))

################################################################################
# For timeseries analysis (pt1):

variety <- as.factor(s2c$variety)
full_design <- model.matrix(~ variety)

################################################################################
# Prepare the sleuth object

so <- sleuth_prep(s2c,
                  full_model = full_design,
                  target_mapping = ttg,
                  read_bootstrap_tpm=TRUE,
                  extra_bootstrap_summary = TRUE,
                  num_cores=1)

so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta = 'varietyamboise', which_model = 'full')
so <- sleuth_fit(so, formula = ~ variety, fit_name = "reduced")
so <- sleuth_lrt(so, 'reduced', 'full')

################################################################################
# Timeseries analysis (pt2) -- analyse the fit to only keep the sufficiently up/downregulated genes:

lrt_results <- sleuth_results(so, 'varietyamboise', test_type = 'wt')
lrt_results <- as.data.frame(lrt_results)

filtered_lrt_results <- lrt_results |>
  dplyr::filter(b < -2 | b > 2)

################################################################################
# Plot timeseries results:

plot_transcript_heatmap(so, head(filtered_lrt_results, n = 300)$target_id, 'tpm')
dev.off()

################################################################################
# Start shiny browser
sleuth_live(so)

################################################################################
# Make plots

fc_results <- subset(lrt_results,select=c('target_id','b','pval'))
colnames(fc_results) <- c('gene','log2fc','pval')

fc_insign <- dplyr::filter(fc_results,pval>0.05)
fc_sign <- dplyr::filter(fc_results, pval < 0.05)
fc_sign_upreg <- dplyr::filter(fc_results,pval<=0.05,log2fc > 0)
fc_sign_upreg$Expression <- 'Upregulated'
fc_sign_downreg <- dplyr::filter(fc_results,pval<=0.05,log2fc < 0)
dplyr::filter(fc_sign_downreg, -log10(pval) >9)
dplyr::filter(fc_sign_downreg, log2fc < -4)
fc_sign_downreg$Expression <- 'Downregulated'

all_data <- rbind(
  mutate(fc_insign, Expression='Insignificant'),
  fc_sign_upreg,fc_sign_downreg
)

colors <- c('Insignificant' = 'grey', 'Upregulated' = '#3EC184', 'Downregulated' = '#C13E7B')
shapes <- c('Insignificant' = 21, 'Upregulated' = 24, 'Downregulated' = 25)

volc_plt <- ggplot(all_data, aes(log2fc, -log10(pval), color = Expression, fill = Expression, shape = Expression)) +
  geom_point(alpha=0.4) +
  scale_color_manual(values = colors, breaks = c('Upregulated', 'Downregulated')) +
  scale_fill_manual(values = colors, breaks = c('Upregulated', 'Downregulated')) +
  scale_shape_manual(values = shapes, breaks = c('Upregulated', 'Downregulated')) +
  scale_y_continuous(limits=c(-0.1, 30), expand=c(0,0), breaks=seq(0,55,5)) +
  scale_x_continuous(limits=c(-10,10), expand=c(0,0), breaks=seq(-10,10,2.5)) +
  labs(color='Expression') +
  geom_hline(yintercept = max(-log10(fc_insign$pval)), alpha=0.8) +
  geom_vline(xintercept = 0, alpha=0.8) +
  theme_minimal() + 
  xlab('log2(FC)') + ylab('-log10(pVal)') +
  theme(panel.grid.minor = element_blank())

volc_plt
ggsave('volcano_plt_amb.pdf',volc_plt)

fc_plt <- ggplot(all_data, aes(gene, abs(log2fc), color = Expression, fill = Expression, shape = Expression)) +
  geom_point(alpha=0.4) +
  scale_color_manual(values = colors, breaks = c('Upregulated', 'Downregulated')) +
  scale_fill_manual(values = colors, breaks = c('Upregulated', 'Downregulated')) +
  scale_shape_manual(values = shapes, breaks = c('Upregulated', 'Downregulated')) +
  scale_y_continuous(limits=c(0,7), expand=c(0,0), breaks=seq(0,10,1)) +
  labs(color='Expression') +
  geom_hline(yintercept = 2, alpha=0.8) +
  theme_classic() + 
  geom_text_repel(max.overlaps=1, data = subset(all_data, abs(log2fc) > 3.5), aes(label = stringr::str_wrap(gene, 2)), size = 3) +
  ylab('log2(FC)') + xlab('Gene') +
  theme(panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(), axis.text.x=element_blank())

fc_plt

ggsave('log2fc_plt_amb.pdf',fc_plt)
