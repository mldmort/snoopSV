<p align="center">
  <img src="https://github.com/mldmort/snoopSV/blob/main/img/snoopSV.png" width="256" align="right"/>
</p>

# snoopSV
This tool processes long-read WGS data (PacBio and ONT), in the context of Structural Variations (SV), Tandem Repeats (TR) and methylation likelihood extraction.

Given a VCF of input SVs, snoopSV determines number of variant-supporting reads and number of reference-supporting reads for each SV and genotypes the calls with genotype qualities.

Given a set of TR regions in bed format, snoopSV determines the base pair deviation of each spanning read and determines genotypes representing base pair deviations in each phase (given the input reads are phased).

Given a set of genomic regions in bed format, snoopSV extracts methylation likelihood information of each CpG site in the regions for each read in a VCF format.

## Genotype Structural Variations
Given an input VCF file containing SVs (INS, DEL, DUP, INV or BND), snoopSV finds the number of variant-supporting reads and reference-supporting reads in each haplotype (given the input reads are phased with HP tags) from an input BAM file for a single sample, and genotypes the variants. It also asigns a genotyping quality (GQ, probablity of the genotype being correct in Phred scale) and a sample quality (SQ, probablity of the genotype being non-reference in Phred scale) to each variant.

Example command line:
```
snoopsv nontr -v <input vcf> -o <output vcf> -s <sample name> -b <bam file>
```
```
usage: snoopsv nontr [-h] -v VCF_IN -o VCF_OUT -s SAMPLE -b BAM [-c CONTIG] [-n N_SECTION] [-i I_SECTION] [--skip-bed SKIP_BED] [--mapping-quality-thr MAPPING_QUALITY_THR]
                     [--buffer-length BUFFER_LENGTH] [--p-err P_ERR] [--len-ratio-tol LEN_RATIO_TOL] [--ins-len-thr INS_LEN_THR] [--del-len-thr DEL_LEN_THR]
                     [--del-recip-overlap-thr DEL_RECIP_OVERLAP_THR] [--ins-recip-overlap-thr INS_RECIP_OVERLAP_THR] [--dup-recip-overlap-thr DUP_RECIP_OVERLAP_THR]
                     [--inv-recip-overlap-thr INV_RECIP_OVERLAP_THR] [--bnd-pos-tol BND_POS_TOL] [--include-svtype [INCLUDE_SVTYPE ...]] [--exclude-svtype [EXCLUDE_SVTYPE ...]]
                     [--exclude-contig [EXCLUDE_CONTIG ...]]

options:
  -h, --help            show this help message and exit
  -v VCF_IN, --vcf_in VCF_IN
                        input VCF file. required.
  -o VCF_OUT, --vcf_out VCF_OUT
                        output VCF file. required.
  -s SAMPLE, --sample SAMPLE
                        sample of interest, if not present in the input VCF, will be added to it. required.
  -b BAM, --bam BAM     mapped Long read bam file associated with the sample of interest. required.
  -c CONTIG, --contig CONTIG
                        contig name (chrom name) to process. If used, the input VCF index should exist. if not present process all contigs.
  -n N_SECTION, --n_section N_SECTION
                        split the input VCF file into this number of sections. default=1
  -i I_SECTION, --i_section I_SECTION
                        which section of the input VCF file to process. used to parallelize runs. default=0
  --skip-bed SKIP_BED   skip analyzing calls intersecting with this bed file
  --mapping-quality-thr MAPPING_QUALITY_THR
                        minimum mapping quality for a read to be considered. default=20
  --buffer-length BUFFER_LENGTH
                        number of base pairs to look for reads around an sv. default=500
  --p-err P_ERR         error probability in sv genotyping model. default=0.01
  --len-ratio-tol LEN_RATIO_TOL
                        tolerance for difference between observed sv length in a read and reported sv length in the input vcf to pass the read as supporting. default=0.2
  --ins-len-thr INS_LEN_THR
                        minimum insertion length in CIGAR to be considered. default=20
  --del-len-thr DEL_LEN_THR
                        minimum deletion length in CIGAR to be considered. default=20
  --del-recip-overlap-thr DEL_RECIP_OVERLAP_THR
                        minimum reciprocal overlap for deletions. default=0.8
  --ins-recip-overlap-thr INS_RECIP_OVERLAP_THR
                        minimum reciprocal overlap for insertions. default=0.5
  --dup-recip-overlap-thr DUP_RECIP_OVERLAP_THR
                        minimum reciprocal overlap for duplications. default=0.8
  --inv-recip-overlap-thr INV_RECIP_OVERLAP_THR
                        minimum reciprocal overlap for inversions. default=0.8
  --bnd-pos-tol BND_POS_TOL
                        tolerance for breakend location in base pairs. default=50
  --include-svtype [INCLUDE_SVTYPE ...]
                        space separated svtypes you want to include
  --exclude-svtype [EXCLUDE_SVTYPE ...]
                        space separated svtypes you want to exclude
  --exclude-contig [EXCLUDE_CONTIG ...]
                        contig names to exclude from processing.
```
## Genotype Tandem Repeats
Given a BED file containing TR regions of interest, snoopSV detects the base pair deviation from the reference for each spanning long read in each TR region and reports them in the output VCF file. This data is reported for each haplotype separately in the VCF file (given the input reads are phased with HP tags). A genotype is assigned based on the median of the base pair deviations in the reads for each haplotype. Currently, there is no genotyping quality reported, however, the number of reads (in each haplotype) is a metric that can be used for filteration.

Example command line:
```
snoopsv tr -a <tr bed file> -v <input vcf> -o <output vcf> -s <sample name> -b <bam file>
```

### Notes about the input parameters
- The purpose of the \<input vcf\> is for its header to be copied to the output file. So you can have your default header lines (such as contig definitions) in this file. Other necessary INFO and FORMAT fields will be added to the output header by snoopSV. There should not be any variants or any samples in this VCF file.
- \<tr bed file\> should have at least three columns: chrom, start, end. If extra columns exist, you should name all columns in the first comment line (starting with "#"), or alternatively provide the column names with an input option: `--annot-columns`. The values of the extra columns will be added to the INFO field of the output VCF unless otherwise is explicitly specified with options such as: `--include-columns` or `--exclude-columns`.

```
usage: snoopsv tr [-h] -a ANNOT_FILE -v VCF_IN -o VCF_OUT [--header-lines HEADER_LINES] -s SAMPLE -b BAM [-c CONTIG] [-n N_SECTION] [-i I_SECTION] [-l L_MAX] [-r R_MIN]
                  [--mapping-quality-thr MAPPING_QUALITY_THR] [--buffer-length BUFFER_LENGTH] [--annot-columns ANNOT_COLUMNS] [--include-columns INCLUDE_COLUMNS]
                  [--exclude-columns EXCLUDE_COLUMNS] [--annot-len-dev-column ANNOT_LEN_DEV_COLUMN] [--annot-len-dev-prc ANNOT_LEN_DEV_PRC] [--annot-len-dev-bp ANNOT_LEN_DEV_BP]
                  [--info-prefix INFO_PREFIX] [--skip-bed SKIP_BED] [--flanking-bp FLANKING_BP]

options:
  -h, --help            show this help message and exit
  -a ANNOT_FILE, --annot_file ANNOT_FILE
                        annotation file. Should be an unzipped bed file with the first three columns: chrom, start, end. If extra columns exist, they will be added to the INFO field of
                        the output VCF. Must have a # header to name the extra fields or alternatively provide the header names with another input option.
  -v VCF_IN, --vcf_in VCF_IN
                        an input VCF file only to use for its header. Make sure no extra samples are present. The presence of the main sample is optional, and will be added if it is not.
  -o VCF_OUT, --vcf_out VCF_OUT
                        output VCF file
  --header-lines HEADER_LINES
                        a file with header lines wants to be added to the output VCF file.
  -s SAMPLE, --sample SAMPLE
                        sample of interest, if not present in the input VCF will be added to it
  -b BAM, --bam BAM     mapped Long read bam file associated with the sample of interest
  -c CONTIG, --contig CONTIG
                        contig name. If used the input VCF should be indexed
  -n N_SECTION, --n_section N_SECTION
                        number of sections in the input VCF file
  -i I_SECTION, --i_section I_SECTION
                        which section of the input VCF file to process
  -l L_MAX, --l_max L_MAX
                        maximum TR span for genotyping
  -r R_MIN, --r_min R_MIN
                        minimum number of reads required in a haplotype for genotyping
  --mapping-quality-thr MAPPING_QUALITY_THR
                        minimum mapping quality for a read to be considered
  --buffer-length BUFFER_LENGTH
                        number of base pairs to look for reads around an annotation
  --annot-columns ANNOT_COLUMNS
                        comma-separated column names of the annotation file.
  --include-columns INCLUDE_COLUMNS
                        comma-separated extra columns you want to include. Default: all
  --exclude-columns EXCLUDE_COLUMNS
                        comma-separated extra columns you want to exclude. Default: None
  --annot-len-dev-column ANNOT_LEN_DEV_COLUMN
                        column name of length deviation or comma-separated length deviations for the annotation. The minimun value of the comma-separated ones is used for significant
                        deviation. Default: None
  --annot-len-dev-prc ANNOT_LEN_DEV_PRC
                        this percentage of annotation length is used as length deviation threshold for significant deviation. Default: None
  --annot-len-dev-bp ANNOT_LEN_DEV_BP
                        number of base pairs to be used as length deviation threshold for significant deviation. Default: None
  --info-prefix INFO_PREFIX
                        prefix of extra columns in the info field. Default: TR
  --skip-bed SKIP_BED   skip analyzing calls intersecting with this bed file
  --flanking-bp FLANKING_BP
                        number of flanking base pairs around an annotation region to match the reads
```
## Extract methylation likelihoods



to further inspect already called SVs (by an SV caller such as sniffles) from long read data, and add more information to the FORMAT feild. More information  includes:
1. Number of reads supporting the variant.
2. Number of reads supporting the reference.
3. Number of reads supporting the variant/reference in each haplotype.
4. Phased/unphased genotype and genotype likelihood.

## What you nead:
1. Structural variant calls in a VCF format.
2. Mapped (optionally haplotagged) bam files for each sample.

## Installation
Install with conda:
```
git clone https://github.com/mldmort/snoopSV.git
cd snoopSV
conda env create -f requirements.txt -n <snoopsv>
conda activate <snoopsv>
make
```
Install with micromamba:
```
git clone https://github.com/mldmort/snoopSV.git
cd snoopSV
micromamba env create -n <snoopsv> -f requirements.txt python=3.10
micromamba activate <snoopsv>
make
```
## How to run snoopSV:
```
snoopsv nontr -v <input vcf> -o <output vcf> -s <sample-bam file>
```
`<sample-bam file>` is a tab delimitted file with the first column being the sample name and the second column the path to the sample mapped bam file.
Each line represents one sample in the input vcf file.
