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
	tr_annot_file = args.tr_annot
	vcf_in = args.vcf_in
	vcf_out = args.vcf_out
	sample_bam_file = args.sample_bam_file
	chrom = args.contig
	n_sec = args.n_section
	i_sec = args.i_section
	tr_span_max = args.l_max
	r_min = args.r_min
	for x, y in args.__dict__.items():
		print(x,':', y)
	sys.stdout.flush()

	GT_TR(tr_annot_file=tr_annot_file, vcf_in=vcf_in, vcf_out=vcf_out, contig=chrom, sample_bam_file=sample_bam_file, n_sec=n_sec, i_sec=i_sec, tr_span_max=tr_span_max, r_min=r_min, verbose=1)

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

	parser_nontr = subparsers.add_parser('nontr', help='Analyze SVs in the non-TR mode')
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

	parser_tr = subparsers.add_parser('tr', help='process TR annotations')
	parser_tr.add_argument('-t', '--tr_annot', required=True, help='TR annotation file. Should be a tab delimited file with columns: tr_chrom, tr_start, tr_end, period length, copy number, period sequence. If you do not have any column information you should put dummy strings for them.')
	parser_tr.add_argument('-v', '--vcf_in', required=True, help='dummy input VCF file to use for header')
	parser_tr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_tr.add_argument('-s', '--sample_bam_file', required=True, help='a map between samples and bam files as a tab delimited text file. First column is the samples, and second column is the absolute path to the bam files')
	parser_tr.add_argument('-c', '--contig', default=None, help='contig name. If used the input VCF should be indexed')
	parser_tr.add_argument('-n', '--n_section', type=int, default=1, help='number of sections in the input VCF file')
	parser_tr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
	parser_tr.add_argument('-l', '--l_max', type=int, default=int(1e9), help='maximum TR span for genotyping')
	parser_tr.add_argument('-r', '--r_min', type=int, default=int(1), help='minimum number of reads required in a haplotype for genotyping')
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
