
def add_header_lines(header_in):

	new_header_INFO = [
	'##INFO=<ID=SKIP_REGION,Number=0,Type=Flag,Description="Skip this call for genotyping, since it falls into skipped regions">',
	'##INFO=<ID=TR_ANNOT,Number=0,Type=Flag,Description="Tandem Repeat annotation genotyped, from genSV">',
	'##INFO=<ID=TR,Number=0,Type=Flag,Description="TR region, target for TR genotyping">',
	'##INFO=<ID=TR_REPEAT_LEN,Number=.,Type=String,Description="TR repeat length">',
	'##INFO=<ID=TR_REPEAT_SEQ,Number=.,Type=String,Description="TR repeat sequence">',
	'##INFO=<ID=TR_REPEAT_START,Number=.,Type=String,Description="TR repeat start">',
	'##INFO=<ID=TR_REPEAT_END,Number=.,Type=String,Description="TR repeat end">',
	'##INFO=<ID=TR_REPEAT_CN,Number=.,Type=String,Description="TR repeat copy number">'
	]
	new_header_FORMAT = [
	'##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, from genSV">',
	'##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, from genSV">',
	'##FORMAT=<ID=RV_P,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, paternal, from genSV">',
	'##FORMAT=<ID=RR_P,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, paternal, from genSV">',
	'##FORMAT=<ID=RV_M,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, maternal, from genSV">',
	'##FORMAT=<ID=RR_M,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, maternal, from genSV">',
	'##FORMAT=<ID=RV_N,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, not phased, from genSV">',
	'##FORMAT=<ID=RR_N,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, not phased, from genSV">',
	'##FORMAT=<ID=GT_SV,Number=1,Type=String,Description="Genotype of the variant, from genSV">',
	'##FORMAT=<ID=GT_SV_PH,Number=1,Type=String,Description="phased genotype of the variant, from genSV">',
	'##FORMAT=<ID=GQ_SV,Number=1,Type=Integer,Description="Phred-scale genotype quality score, from genSV">',
	'##FORMAT=<ID=SQ_SV,Number=1,Type=Integer,Description="Phred-scale sample quality score, from genSV">',
	'##FORMAT=<ID=P_11,Number=1,Type=Float,Description="probability of 1/1 genotype, from genSV">',
	'##FORMAT=<ID=P_01,Number=1,Type=Float,Description="probability of 0/1 genotype, from genSV">',
	'##FORMAT=<ID=P_00,Number=1,Type=Float,Description="probability of 0/0 genotype, from genSV">',
	'##FORMAT=<ID=GT_TR_AL,Number=.,Type=String,Description="Genotype in TR region, from genSV">',
	'##FORMAT=<ID=GQ_TR_AL,Number=.,Type=String,Description="Phred-scale genotype quality score in TR region, from genSV">',
	'##FORMAT=<ID=CN_TR_AL_H1,Number=.,Type=String,Description="Counts of repeat unit in TR region for H1, from genSV">',
	'##FORMAT=<ID=CN_TR_AL_H2,Number=.,Type=String,Description="Counts of repeat unit in TR region for H2, from genSV">',
	'##FORMAT=<ID=CN_TR_AL_H0,Number=.,Type=String,Description="Counts of repeat unit in TR region for H0, from genSV">',
	'##FORMAT=<ID=GT_TR_LN,Number=.,Type=String,Description="Genotype in TR region, from genSV">',
	'##FORMAT=<ID=GQ_TR_LN,Number=.,Type=String,Description="Phred-scale genotype quality score in TR region, from genSV">',
	'##FORMAT=<ID=CN_TR_LN_H1,Number=.,Type=String,Description="Counts of repeat unit in TR region for H1, from genSV">',
	'##FORMAT=<ID=CN_TR_LN_H2,Number=.,Type=String,Description="Counts of repeat unit in TR region for H2, from genSV">',
	'##FORMAT=<ID=CN_TR_LN_H0,Number=.,Type=String,Description="Counts of repeat unit in TR region for H0, from genSV">',
	'##FORMAT=<ID=GT_TR_BP,Number=.,Type=String,Description="Genotype in TR region, from genSV">',
	'##FORMAT=<ID=BP_DV,Number=1,Type=String,Description="base pair deviation from reference in TR region, from genSV">',
	'##FORMAT=<ID=CN_TR_BP_H1,Number=.,Type=String,Description="Counts of base pairs in TR region for H1, from genSV">',
	'##FORMAT=<ID=CN_TR_BP_H2,Number=.,Type=String,Description="Counts of base pairs in TR region for H2, from genSV">',
	'##FORMAT=<ID=CN_TR_BP_H0,Number=.,Type=String,Description="Counts of base pairs in TR region for H0, from genSV">'
	]

	for line in new_header_INFO + new_header_FORMAT:
		header_in.add_line(line)

	return
