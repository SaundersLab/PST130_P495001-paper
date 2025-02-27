# Scripts used for data analysis in the _PST130_P495001_ paper
## Count SNPs per base
After running variant calling and creating consensus sequences of the transcriptomic samples ([STAR](https://github.com/alexdobin/STAR) version 2.5; [Samtools](https://www.htslib.org/) version 0.1.19), using _Pst_ isolate 104E137A-[^1], use the `snpb_script.py` to extract the average and standard deviation SNP/b, for the housekeeping genes and _PST130_P495001_.

## Identify Variants/Haplotypes/Isoforms
After running variant calling and creating consensus sequences of the transcriptomic samples ([STAR](https://github.com/alexdobin/STAR) version 2.5; [Samtools](https://www.htslib.org/) version 0.1.19), using _Pst_ isolate 104E137A-[^1], use the `gene_variant_identify.ipynb` to:
+ Find the different haplotype variants of the gene of interest, given no ambiguities are present in the consensus sequence.
+ Find the different isoforms from the gene of interest.
+ Align the different haplotypes and save a `.png` image of the MSA alignment.

## Extract information from the Expression Browser[^2]
To identify which isolates of the ones included in the expression browser, lack expression of _PST130_P495001_, the data was downloaded (directory: REB_for_figshare) and everything apart from the Dobon et al. samples (as they are from a timeseries) were used in the analysis. The transcript counts were used to remove isolates that had a low expression of the housekeeping genes, as a way of filtering to include only high quality reads. For this, any isolates with a mean expression <15% the median were removed. Out of the 888 remaining isolates, 10 were identified to lack expression of _PST130_P495001_. This data analysis was carried out using `check_expression.ipynb`.

## Running Differential Gene Expression Analysis
To obtain the abundance counts of the transcriptomic reads, pseudoalignments were performed using [Kallisto](https://pachterlab.github.io/kallisto/) version 0.51.1. The transcript counts were then analysed using [Sleuth](https://pachterlab.github.io/sleuth/) version 0.30.1, using `sleuth.R`.


[^1]: Schwessinger, B. _et al._ A Near-Complete Haplotype-Phased Genome of the Dikaryotic Wheat Stripe Rust Fungus Puccinia striiformis f. sp. tritici Reveals High Interhaplotype Diversity. _mBio_ **9** (2018). [10.1128/mBio.02275-17](https://doi.org:10.1128/mBio.02275-17)
[^2]: [rust-expression.com](https://www.rust-expression.com/)
