<p align="center">
  <img src="https://github.com/mldmort/snoopSV/blob/main/img/snoopSV.png" width="256" align="right"/>
</p>

# snoopSV
This tool can be used to further inspect already called SVs (by an SV caller such as sniffles) from long read data, and add more information to the FORMAT feild. More information  includes:
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
micromamba env create -n <snoopsv> python=3.10
micromamba activate <snoopsv>
make
```
## How to run snoopSV:
```
snoopsv nontr -v <input vcf> -o <output vcf> -s <sample-bam file>
```
`<sample-bam file>` is a tab delimitted file with the first column being the sample name and the second column the path to the sample mapped bam file.
Each line represents one sample in the input vcf file.
