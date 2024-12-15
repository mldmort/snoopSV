import argparse
import sys
from snoopsv.snoop_sv import GT_nonTR
from snoopsv.snoop_tr import GT_TR
from snoopsv.snoop_score import SCORE_VCF

def run_nontr(args):
	vcf_in = args.vcf_in
	vcf_out = args.vcf_out
	sample = args.sample
	bam = args.bam
	chrom = args.contig
	n_sec = args.n_section
	i_sec = args.i_section
	skip_bed = args.skip_bed
	mapping_quality_thr = args.mapping_quality_thr
	buffer_length = args.buffer_length
	p_err = args.p_err
	len_ratio_tol = args.len_ratio_tol
	ins_len_thr = args.ins_len_thr
	del_len_thr = args.del_len_thr
	del_recip_overlap_thr = args.del_recip_overlap_thr
	ins_recip_overlap_thr = args.ins_recip_overlap_thr
	dup_recip_overlap_thr = args.dup_recip_overlap_thr
	inv_recip_overlap_thr = args.inv_recip_overlap_thr
	bnd_pos_tol = args.bnd_pos_tol
	include_svtype = args.include_svtype
	exclude_svtype = args.exclude_svtype
	exclude_contig = args.exclude_contig
	for x, y in args.__dict__.items():
		print(x,':', y)
	sys.stdout.flush()

	GT_nonTR(vcf_in, vcf_out, contig=chrom, sample=sample, bam=bam, n_sec=n_sec, i_sec=i_sec,
			 skip_bed=skip_bed, mapping_quality_thr=mapping_quality_thr, buffer_length=buffer_length,
			 p_err=p_err, len_ratio_tol=len_ratio_tol, ins_len_thr=ins_len_thr, del_len_thr=del_len_thr,
			 del_recip_overlap_thr=del_recip_overlap_thr, ins_recip_overlap_thr=ins_recip_overlap_thr,
			 dup_recip_overlap_thr=dup_recip_overlap_thr, inv_recip_overlap_thr=inv_recip_overlap_thr,
			 bnd_pos_tol=bnd_pos_tol, verbose=1, include_svtype=include_svtype, exclude_svtype=exclude_svtype,
			 exclude_contig=exclude_contig)

def run_tr(args):
	annot_file = args.annot_file
	vcf_in = args.vcf_in
	vcf_out = args.vcf_out
	sample = args.sample
	bam = args.bam
	chrom = args.contig
	n_sec = args.n_section
	i_sec = args.i_section
	tr_span_max = args.l_max
	r_min = args.r_min
	mapping_quality_thr = args.mapping_quality_thr
	buffer_length = args.buffer_length
	annot_columns = args.annot_columns
	include_columns = args.include_columns
	exclude_columns = args.exclude_columns
	annot_len_dev_column = args.annot_len_dev_column
	annot_len_dev_prc = args.annot_len_dev_prc
	annot_len_dev_bp = args.annot_len_dev_bp
	header_lines = args.header_lines
	info_prefix = args.info_prefix
	skip_bed = args.skip_bed
	flanking_bp = args.flanking_bp
	for x, y in args.__dict__.items():
		print(x,':', y)
	sys.stdout.flush()

	GT_TR(annot_file=annot_file, vcf_in=vcf_in, vcf_out=vcf_out, contig=chrom, sample=sample, bam=bam, n_sec=n_sec, i_sec=i_sec, mapping_quality_thr=mapping_quality_thr, buffer_length=buffer_length, tr_span_max=tr_span_max, r_min=r_min, annot_columns=annot_columns, include_columns=include_columns, exclude_columns=exclude_columns, annot_len_dev_column=annot_len_dev_column, annot_len_dev_prc=annot_len_dev_prc, annot_len_dev_bp=annot_len_dev_bp, header_lines=header_lines, info_prefix=info_prefix, skip_bed=skip_bed, flanking_bp=flanking_bp, verbose=1)

def run_score(args):
	vcf_in = args.vcf_in
	annot_in = args.annot_in
	cov_in = args.coverage_file
	vcf_out = args.vcf_out
	models = args.models
	for x, y in args.__dict__.items():
		print(x,':', y)

	SCORE_VCF(vcf_in=vcf_in, annot_in=annot_in, cov_in=cov_in, models=models, vcf_out=vcf_out)

def main():
	parser = argparse.ArgumentParser(description='snoopSV detects long reads supporting SV calls, and genotype them.')
	subparsers = parser.add_subparsers(help='available sub-commands')

	parser_nontr = subparsers.add_parser('nontr', help='Detect SV signitures in non-TR mode')
	parser_nontr.add_argument('-v', '--vcf_in', required=True, help='input VCF file')
	parser_nontr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_nontr.add_argument('-s', '--sample', required=True, help='sample of interest, if not present in the input VCF will be added to it')
	parser_nontr.add_argument('-b', '--bam', required=True, help='mapped Long read bam file associated with the sample of interest')
	parser_nontr.add_argument('-c', '--contig', default=None, help='contig name. If used, the input VCF index should exist')
	parser_nontr.add_argument('-n', '--n_section', type=int, default=1, help='split the input VCF file into this number of sections')
	parser_nontr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
	parser_nontr.add_argument('--skip-bed', default=None, help='skip analyzing calls intersecting with this bed file')
	parser_nontr.add_argument('--mapping-quality-thr', default=20, type=int, help='minimum mapping quality for a read to be considered')
	parser_nontr.add_argument('--buffer-length', default=500, type=int, help='number of base pairs to look for reads around an sv')
	parser_nontr.add_argument('--p-err', default=0.01, type=float, help='error probability in sv genotyping model')
	parser_nontr.add_argument('--len-ratio-tol', default=0.2, type=float, help='tolerance for difference between observed sv length in a read and reported sv length in the input vcf to pass the read as supporting')
	parser_nontr.add_argument('--ins-len-thr', default=20, type=int, help='minimum insertion length in CIGAR to be considered')
	parser_nontr.add_argument('--del-len-thr', default=20, type=int, help='minimum deletion length in CIGAR to be considered')
	parser_nontr.add_argument('--del-recip-overlap-thr', default=0.8, type=float, help='minimum reciprocal overlap for deletions')
	parser_nontr.add_argument('--ins-recip-overlap-thr', default=0.5, type=float, help='minimum reciprocal overlap for insertions')
	parser_nontr.add_argument('--dup-recip-overlap-thr', default=0.8, type=float, help='minimum reciprocal overlap for duplications')
	parser_nontr.add_argument('--inv-recip-overlap-thr', default=0.8, type=float, help='minimum reciprocal overlap for inversions')
	parser_nontr.add_argument('--bnd-pos-tol', default=50, type=int, help='tolerance for breakend location in base pairs')
	parser_nontr.add_argument('--include-svtype', default=None, type=str, nargs='*', help='space separated svtypes you want to include')
	parser_nontr.add_argument('--exclude-svtype', default=None, type=str, nargs='*', help='space separated svtypes you want to exclude')
	parser_nontr.add_argument('--exclude-contig', default=None, type=str, nargs='*', help='contig names to exclude from processing.')
	parser_nontr.set_defaults(func=run_nontr)

	parser_tr = subparsers.add_parser('tr', help='Detect sequence expansion/contraction across TR annotations')
	parser_tr.add_argument('-a', '--annot_file', required=True, help='annotation file. Should be an unzipped bed file with the first three columns: chrom, start, end. If extra columns exist, they will be added to the INFO field of the output VCF. Must have a # header to name the extra fields or alternatively provide the header names with another input option.')
	parser_tr.add_argument('-v', '--vcf_in', required=True, help='an input VCF file only to use for its header. Make sure no extra samples are present. The presence of the main sample is optional, and will be added if it is not.')
	parser_tr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_tr.add_argument('--header-lines', default=None, type=str, required=False, help='a file with header lines wants to be added to the output VCF file.')
	parser_tr.add_argument('-s', '--sample', required=True, help='sample of interest, if not present in the input VCF will be added to it')
	parser_tr.add_argument('-b', '--bam', required=True, help='mapped Long read bam file associated with the sample of interest')
	parser_tr.add_argument('-c', '--contig', default=None, help='contig name. If used the input VCF should be indexed')
	parser_tr.add_argument('-n', '--n_section', type=int, default=1, help='number of sections in the input VCF file')
	parser_tr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
	parser_tr.add_argument('-l', '--l_max', type=int, default=int(1e9), help='maximum TR span for genotyping')
	parser_tr.add_argument('-r', '--r_min', type=int, default=int(1), help='minimum number of reads required in a haplotype for genotyping')
	parser_tr.add_argument('--mapping-quality-thr', default=20, type=int, help='minimum mapping quality for a read to be considered')
	parser_tr.add_argument('--buffer-length', default=500, type=int, help='number of base pairs to look for reads around an annotation')
	parser_tr.add_argument('--annot-columns', default=None, type=str, help='comma-separated column names of the annotation file.')
	parser_tr.add_argument('--include-columns', default='all', type=str, help='comma-separated extra columns you want to include. Default: all')
	parser_tr.add_argument('--exclude-columns', default=None, type=str, help='comma-separated extra columns you want to exclude. Default: None')
	parser_tr.add_argument('--annot-len-dev-column', default=None, type=str, help='column name of length deviation or comma-separated length deviations for the annotation. The minimun value of the comma-separated ones is used for significant deviation. Default: None')
	parser_tr.add_argument('--annot-len-dev-prc', default=None, type=float, help='this percentage of annotation length is used as length deviation threshold for significant deviation. Default: None')
	parser_tr.add_argument('--annot-len-dev-bp', default=None, type=int, help='number of base pairs to be used as length deviation threshold for significant deviation. Default: None')
	parser_tr.add_argument('--info-prefix', default='TR', type=str, help='prefix of extra columns in the info field. Default: TR')
	parser_tr.add_argument('--skip-bed', default=None, help='skip analyzing calls intersecting with this bed file')
	parser_tr.add_argument('--flanking-bp', type=int, default=50, help='number of flanking base pairs around an annotation region to match the reads')
	parser_tr.set_defaults(func=run_tr)

	parser_score = subparsers.add_parser('score', help='score variants with a ML model')
	parser_score.add_argument('-v', '--vcf_in', required=True, help='input VCF file, unscored')
	parser_score.add_argument('-a', '--annot_in', required=True, help='input annotation file coresponding to input VCF file')
	parser_score.add_argument('-o', '--vcf_out', required=True, help='output VCF file, scored')
	parser_score.add_argument('-m', '--models', required=True, help='a file with ML model labels in the first column, pickled file addresses in the second column, and camma separated model variables in the third column. each line represent one model to be applied.')
	parser_score.add_argument('-c', '--coverage_file', required=True, help='coverage file, column one: sample name, column two: coverage summary file')
	parser_score.set_defaults(func=run_score)
	### parse the command line
	args = parser.parse_args()

	### run the program
	args.func(args)

if __name__ == '__main__':
	main()
