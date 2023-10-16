import argparse
from snoopsv.snoop_sv import GT_nonTR
from snoopsv.snoop_tr import GT_TR
from snoopsv.snoop_score import SCORE_VCF

def run_nontr(args):
	tr_annot_file = args.tr_annot
	vcf_in = args.vcf_in
	vcf_out = args.vcf_out
	sample_bam_file = args.sample_bam_file
	chrom = args.contig
	n_sec = args.n_section
	i_sec = args.i_section
	for x, y in args.__dict__.items():
		print(x,':', y)
	if tr_annot_file != '':
		print('WARN: tr annotation file provided as an input is NOT utilized in the nonTR genotyper for this version. TR annotation should be done after genotyping.')

	GT_nonTR(vcf_in, vcf_out, contig=chrom, sample_bam_file=sample_bam_file, n_sec=n_sec, i_sec=i_sec, verbose=1)

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
	parser = argparse.ArgumentParser(description='genSV is a genotyper for cohort SV calls from third generation data (ONT/PacBio)')
	subparsers = parser.add_subparsers(help='available sub-commands')

	parser_nontr = subparsers.add_parser('nontr', help='process non-TR SVs')
	parser_nontr.add_argument('-v', '--vcf_in', required=True, help='input VCF file')
	parser_nontr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_nontr.add_argument('-s', '--sample_bam_file', required=True, help='a map between samples and bam files as a tab delimited text file. First column is the samples, and second column is the absolute path to the bam files')
	parser_nontr.add_argument('-t', '--tr_annot', required=False, help='TR annotation file. Should be a tab delimited file with columns: tr_chrom, tr_start, tr_end, period length, copy number, period sequence. If you do not have any column information you should put dummy strings for them.', default='')
	parser_nontr.add_argument('-c', '--contig', default=None, help='contig name. If used the input VCF should be indexed')
	parser_nontr.add_argument('-n', '--n_section', type=int, default=1, help='number of sections in the input VCF file')
	parser_nontr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
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
