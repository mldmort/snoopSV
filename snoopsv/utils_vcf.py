
def add_header_lines(header_in):

	new_header_INFO = [
	'##INFO=<ID=SKIP_REGION,Number=0,Type=Flag,Description="Skip this call for genotyping">'
	]
	new_header_FORMAT = [
	'##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, from snoopSV">',
	'##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, from snoopSV">',
	'##FORMAT=<ID=RV_P,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, paternal, from snoopSV">',
	'##FORMAT=<ID=RR_P,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, paternal, from snoopSV">',
	'##FORMAT=<ID=RV_M,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, maternal, from snoopSV">',
	'##FORMAT=<ID=RR_M,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, maternal, from snoopSV">',
	'##FORMAT=<ID=RV_N,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, not phased, from snoopSV">',
	'##FORMAT=<ID=RR_N,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, not phased, from snoopSV">',
	'##FORMAT=<ID=GT_SV,Number=1,Type=String,Description="Genotype of the variant, from snoopSV">',
	'##FORMAT=<ID=GT_SV_PH,Number=1,Type=String,Description="phased genotype of the variant, from snoopSV">',
	'##FORMAT=<ID=GQ_SV,Number=1,Type=Integer,Description="Phred-scale genotype quality score, from snoopSV">',
	'##FORMAT=<ID=SQ_SV,Number=1,Type=Integer,Description="Phred-scale sample quality score, from snoopSV">',
	'##FORMAT=<ID=P_11,Number=1,Type=Float,Description="probability of 1/1 genotype, from snoopSV">',
	'##FORMAT=<ID=P_01,Number=1,Type=Float,Description="probability of 0/1 genotype, from snoopSV">',
	'##FORMAT=<ID=P_00,Number=1,Type=Float,Description="probability of 0/0 genotype, from snoopSV">'
	]

	for line in new_header_INFO + new_header_FORMAT:
		header_in.add_line(line)

	return

def add_header_lines_tr(header_in, prefix='', extra_fields_list=None):

	new_header_INFO = [
	f'##INFO=<ID=SKIP_REGION,Number=0,Type=Flag,Description="Skip this call for genotyping">',
	f'##INFO=<ID={prefix}_ANNOT,Number=0,Type=Flag,Description="{prefix} annotation genotyped, from snoopSV">',
	]
	if extra_fields_list:
		for extra_field in extra_fields_list:
			line = f'##INFO=<ID={prefix}_{extra_field},Number=.,Type=String,Description="{extra_field} from annotation file, from snoopSV">'
			new_header_INFO.append(line)
	new_header_FORMAT = [
	'##FORMAT=<ID=GT_BP,Number=1,Type=String,Description="phased number of base pairs across annotations, from snoopSV">',
	'##FORMAT=<ID=GT_DEV,Number=1,Type=String,Description="phased number of bp deviation across annotations, from snoopSV">',
	'##FORMAT=<ID=BP_H1,Number=.,Type=String,Description="number of base pairs for each read in H1, from snoopSV">',
	'##FORMAT=<ID=BP_H2,Number=.,Type=String,Description="number of base pairs for each read in H2, from snoopSV">',
	'##FORMAT=<ID=BP_H0,Number=.,Type=String,Description="number of base pairs for each read in H0, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H1,Number=.,Type=String,Description="number of bp deviation for each read in H1, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H2,Number=.,Type=String,Description="number of bp deviation for each read in H2, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H0,Number=.,Type=String,Description="number of bp deviation for each read in H0, from snoopSV">',
	'##FORMAT=<ID=N_H1,Number=1,Type=Integer,Description="number of supporting reads in H1, from snoopSV">',
	'##FORMAT=<ID=N_H2,Number=1,Type=Integer,Description="number of supporting reads in H2, from snoopSV">',
	'##FORMAT=<ID=N_H0,Number=1,Type=Integer,Description="number of supporting reads in H0, from snoopSV">',
	'##FORMAT=<ID=SA_SUPP,Number=.,Type=String,Description="in tr mode, if any read signature comes from supplementary alignments, 0: no, 1: yes">',
	#'##FORMAT=<ID=GQ_H1,Number=.,Type=String,Description="quality of GT_BP in H1, from snoopSV">',
	#'##FORMAT=<ID=GQ_H2,Number=.,Type=String,Description="quality of GT_BP in H2, from snoopSV">',
	#'##FORMAT=<ID=GQ_H0,Number=.,Type=String,Description="quality of GT_BP in H0, from snoopSV">',
	]
	for line in new_header_INFO + new_header_FORMAT:
		header_in.add_line(line)

	return

def add_header_lines_methyl(header_in, prefix='', extra_fields_list=None):

	new_header_INFO = [
	f'##INFO=<ID=SKIP_REGION,Number=0,Type=Flag,Description="Skip this call for genotyping">',
	]
	if extra_fields_list:
		for extra_field in extra_fields_list:
			line = f'##INFO=<ID={prefix}_{extra_field},Number=.,Type=String,Description="{extra_field} from annotation file, from snoopSV">'
			new_header_INFO.append(line)
	new_header_FORMAT = [
	'##FORMAT=<ID=METHYL_H1,Number=.,Type=String,Description="methylation probabilities, 0-100, for each Cm reported for each read in H1, from snoopSV">',
	'##FORMAT=<ID=METHYL_H2,Number=.,Type=String,Description="methylation probabilities, 0-100, for each Cm reported for each read in H2, from snoopSV">',
	'##FORMAT=<ID=METHYL_H0,Number=.,Type=String,Description="methylation probabilities, 0-100, for each Cm reported for each read in H0, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H1,Number=.,Type=String,Description="number of bp deviation for each read in H1, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H2,Number=.,Type=String,Description="number of bp deviation for each read in H2, from snoopSV">',
	'##FORMAT=<ID=BP_DEV_H0,Number=.,Type=String,Description="number of bp deviation for each read in H0, from snoopSV">',
	'##FORMAT=<ID=RN_H1,Number=.,Type=String,Description="read names in H1, from snoopSV">',
	'##FORMAT=<ID=RN_H2,Number=.,Type=String,Description="read names in H2, from snoopSV">',
	'##FORMAT=<ID=RN_H0,Number=.,Type=String,Description="read names in H0, from snoopSV">',
	'##FORMAT=<ID=N_H1,Number=1,Type=Integer,Description="number of supporting reads in H1, from snoopSV">',
	'##FORMAT=<ID=N_H2,Number=1,Type=Integer,Description="number of supporting reads in H2, from snoopSV">',
	'##FORMAT=<ID=N_H0,Number=1,Type=Integer,Description="number of supporting reads in H0, from snoopSV">',
	'##FORMAT=<ID=PS_H1,Number=.,Type=String,Description="phase set for read names in H1, from snoopSV">',
	'##FORMAT=<ID=PS_H2,Number=.,Type=String,Description="phase set for read names in H2, from snoopSV">',
	]
	for line in new_header_INFO + new_header_FORMAT:
		header_in.add_line(line)

	return
