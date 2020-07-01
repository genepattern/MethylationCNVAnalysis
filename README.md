# MethylationCNVAnalysis v2.4
## Genome-wide copy number estimation and segmentation from Illumina 450k or EPIC Methylation arrays

This is a GenePattern module written in R v3.6.1.

Authors: Clarence Mah, Owen Chapman


Contact:

https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!forum/genepattern-help

ochapman@ucsd.edu


Algorithm Version: Conumee 1.0.0

## Summary
The Illumina Methylation Array assay family queries methylation levels at 450,000 (450k) or 850,000 (EPIC) locations across the human genome. This module estimates copy number variation (CNV) of one or more tumor samples by comparing methylation levels at these locations to a reference set of "normal" methylation profiles.

### Requires:
- Tumor methylation profiles: [sample]_Grn.idat and [sample]_Red.idat
- A set of control methylation profiles in the same format.

## References
Mah CK, Mesirov JP, Chavez L. An accessible GenePattern notebook for the copy number variation analysis of Illumina Infinium DNA methylation arrays. F1000Res. 2018;7:ISCB Comm J-1897. Published 2018 Dec 5. doi:10.12688/f1000research.16338.1

Aryee MJ, Jaffe AE, Corrada-Bravo H, et al. Minfi: a flexible and comprehensive Bioconductor package for the analysis of Infinium DNA methylation microarrays. Bioinformatics. 2014;30(10):1363–1369. doi:10.1093/bioinformatics/btu049

Hovestadt V, Zapatka M. conumee: Enhanced copy-number variation analysis using Illumina DNA methylation arrays. R package version 1.9.0, http://bioconductor.org/packages/conumee/.

Seshan VE, Olshen A (2019). DNAcopy: DNA copy number data analysis. R package version 1.58.0.

## Parameters
| Name  |	Description |
|-------|-------------|
| query sample data *	| A ZIP or GZ file containing your sample(s) methylation microarray data in the Illumina Demo Dataset folder structure.<br>This folder should contain matched green and red channel files for each query sample, formatted as \[sample\]\_Grn.idat and \[sample\]\_Red.idat ; then compressed into .tar.gz or .zip. |
| control sample names | Specify a list of control sample names separated by commas if the control sample data is included in the "sample data" file. Otherwise write "none" if control samples data is separate. |
| control sample data *	| A ZIP or GZ file containing your sample(s) methylation microarray data in the Illumina Demo Dataset folder structure. This folder should contain matched green and red channel files for each query sample, formatted as \[sample\]\_Grn.idat and \[sample\]\_Red.idat ; then compressed into .tar.gz or .zip.<br>An example control set, consisting of 72 noncancer central nervous system samples, may be found at https://datasets.genepattern.org/data/module_support_files/MethylationCNVAnalysis/Methyl_cnv/CNS_450k_controls.tar.gz. Appropriate control sets may also be obtained from the ENCODE database. |
| genes to highlight	| A file with a list of genes (HUGO gene symbols) to highlight in output plots. Format as one gene symbol per line. An example file of common cancer genes may be found at [ftp://gpftp.broadinstitute.org/methylation/common_cancer_genes.txt](ftp://gpftp.broadinstitute.org/methylation/common_cancer_genes.txt) |
| ignore regions	| .Bed file: regions which should be excluded from CN analysis (hg19). A standard file of highly repetitive regions may be found at [ftp://gpftp.broadinstitute.org/methylation/ignore_regions.bed](ftp://gpftp.broadinstitute.org/methylation/ignore_regions.bed) |
| sex chromosomes *	| Include CN estimates of X & Y chromosomes. |

\* required

## Output Files
Name |	Description
-----|------------
[sample].cnv.seg	| Tab-separated segmentation file of copy number segments and their estimated log2 fold change from controls, output by DNACopy. Format: ID    chrom    loc.start    loc.end    num.mark    bstat    pval    seg.mean    seg.median
[sample].cnvPlots.pdf	| Copy number plots output by conumee. Contains genome-wide copy number estimations, as well as plots at each gene of interest.
[sample].detail.cnv.seg	| Tab-separated file of copy number estimates at each gene of interest.<br>Format: chr    start    end    name    sample    probes.Var1    probes.Freq    value
qcPlots.pdf	| QC plots output by conumee. Check to ensure that your controls were appropriate.

## License
MethylationCNVAnalysis is distributed under a modified BSD license available at https://raw.githubusercontent.com/genepattern/methylation_cnv_analysis_module/master/LICENSE

### Platform Dependencies
Task Type:


CPU Type:
any

Operating System:
any

Language:
R v3.4.1

Versions:
2.4 | 06-2020 | Geneslist of detail regions now an optional argument
----|---------|----------------
2.3 | 11-2019 | Controls can be EPIC 2019 arrays
2.2 | 09-2019 | Minor text fixes
2.1 | 09-2019 | Minor text fixes
2.0 | 09-2019 | EPIC 2019 update
1.0 | 08-2018 | Initial release
