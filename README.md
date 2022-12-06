# Some scripts for genomics (_in construction_)

## Geral information about the scripts:

###### File [_script_bwa_UCEs_2022.sh_](https://github.com/SergioB1983/some_scripts_for_genomics/blob/main/script_bwa_UCEs_2022.sh):

In this file, there are a very simple series of commands that allow to create `.bam` files needed to implement the phasing of UCEs data, that includes the extraction and construction of a SNPs (Single Nucleotide polymorphisms) matrix. This methodology was proposed by Harvey et al. (2016) and includes some modifications by Glaucia Del-Rio (Cornell Laboratory of Ornithology). 

The commands use **samtools**, **bwa-mem**, and **gatk** scripts, and does the same processes included in the script `phyluce_snp_bwa_multiple_align` from the software phyluce (v.1.6). This script is not included in the latest version (v.1.7). 

For more information, please read the script!

**Bibliography**

Faircloth BC. 2015. PHYLUCE is a software package for the analysis of conserved genomic loci. Bioinformatics. doi: 10.1093/bioinformatics/btv646.

Harvey, M. G., Smith, B. T., Glenn, T. C., Faircloth, B. C., Brumfield, R. T. (2016). Sequence Capture versus Restriction Site Associated DNA Sequencing for Shallow Systematics. Systematic Biology, 65(5), 910–924. https://doi.org/10.1093/sysbio/syw036.

