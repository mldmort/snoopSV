import pysam
import glob
import pandas as pd
import numpy as np
import sys
from math import log10, factorial
import subprocess
from genSV import sv_class, infer_gt, infer_gt_str, sv_signiture, str_signature

def make_VCF_GT(vcf_in, vcf_out, contig, bam_file, n_sec, i_sec, verbose=1):

	### genotyping setting
	mapping_quality_thr = 20
	region_buffer_length = 1000
	p_err = 0.05

	STR_file_TRF_target = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot_per-2-250.bed'
	STR_file_TRF_all = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed'
	STR_file_RM_Simpe = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Repeats_Masker_Simple_repeat.bed'
	STR_file_MG = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/hg38_ver13.bed'

	REF_FILE = '/oasis/scratch/comet/smmortazavi/temp_project/HUMAN_DATA/REF_GENOME/GRCh38_full_analysis_set_plus_decoy_hla.fa'
	fa_handle = pysam.FastaFile(REF_FILE)

	fh_bam = pysam.AlignmentFile(bam_file, 'rb')

	k_s_dict = {}
	p_val_thr = 0.07
	for k in range(2,262):
		for s_sel in range(k,0,-1):
			p_val = factorial(k)/factorial(s_sel)/factorial(k-s_sel)/float(4**s_sel)
			if p_val > p_val_thr:
				k_s_dict[k] = s_sel + 1
				break

	if contig:
		command = ('bcftools query -r '+contig+' -f %CHROM\\n '+vcf_in).split()
	else:
		command = ('bcftools query -f %CHROM\\n '+vcf_in).split()
	ps = subprocess.Popen(command, stdout=subprocess.PIPE)
	n_calls = int(subprocess.run(['wc', '-l'], stdin=ps.stdout, check=True, capture_output=True, text=True).stdout)
	### i_sec should go from 0 to n_sec to cover all calls
	i_rec_start = i_sec * int(n_calls / n_sec)
	i_rec_end = (i_sec + 1) * int(n_calls / n_sec)
	if verbose == 1:
		print('n_calls:',  n_calls)
		print('i_sec:', i_sec)
		print('n_sec:', n_sec)
		print('i_rec_start:', i_rec_start)
		print('i_rec_end:', i_rec_end)

	new_header_INFO = ['##INFO=<ID=STR,Number=0,Type=Flag,Description="STR region">',
	'##INFO=<ID=STR_TRF_OTHER,Number=0,Type=Flag,Description="STR region">',
	'##INFO=<ID=STR_RM_SR,Number=0,Type=Flag,Description="STR region">',
	'##INFO=<ID=STR_MG,Number=0,Type=Flag,Description="STR region">',
	'##INFO=<ID=STR_REPEAT_LEN,Number=.,Type=String,Description="STR repeat length, can have overlaping STRs">',
	'##INFO=<ID=STR_REPEAT_SEQ,Number=.,Type=String,Description="STR repeat sequence, can have overlaping STRs">'
	]
	new_header_FORMAT = ['##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, from genSV">',
	'##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, from genSV">',
	'##FORMAT=<ID=GT_SV,Number=1,Type=String,Description="Genotype of the variant, from genSV">',
	'##FORMAT=<ID=GQ_SV,Number=1,Type=Integer,Description="Phred-scale genotype quality score, from genSV">',
	'##FORMAT=<ID=P_11,Number=1,Type=Float,Description="probability of 1/1 genotype, from genSV">',
	'##FORMAT=<ID=P_01,Number=1,Type=Float,Description="probability of 0/1 genotype, from genSV">',
	'##FORMAT=<ID=P_00,Number=1,Type=Float,Description="probability of 0/0 genotype, from genSV">',
	'##FORMAT=<ID=GT_STR,Number=.,Type=String,Description="Genotype in STR region, from genSV. Shows the number of repeats added or deleted">',
	'##FORMAT=<ID=GQ_STR,Number=.,Type=String,Description="Phred-scale genotype quality score in STR region, from genSV">'
	]

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	for line in new_header_INFO+new_header_FORMAT:
		header_in.add_line(line)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	count_skip_region = 0
	count_skip_sec = 0
	for i_rec, rec in enumerate(fh_vcf_in.fetch(contig=contig)):
		if (i_rec < i_rec_start) or (i_rec >= i_rec_end):
			count_skip_sec += 1
			continue
		sv_id = rec.id
		svtype = rec.info['SVTYPE']
		svlen = rec.info['SVLEN']
		chr2 = rec.info['CHR2']
		chrom = rec.chrom
		start = rec.start
		stop = rec.stop

		chr1_region_1=143150000
		chr1_region_2=149900000
		if (chrom=='chr1') and ((start > chr1_region_1) and (start < chr1_region_2)) and ((stop > chr1_region_1) and (stop < chr1_region_2)):
			count_skip_region += 1
			continue
		#print('svtype:', svtype)
		#print('svlen:', svlen) 
		#print('chr2:', chr2) 
		#print('chrom:', chrom)
		#print('start:', start)
		#print('stop:', stop)
		#print('sv_id:', sv_id)
		#print('rec.samples:', rec.samples)
		target_sv = sv_class(svtype, chrom, start, stop, chr2, svlen)

		if svtype == 'INS':
			pos_start = target_sv.start - 50
			pos_stop = target_sv.stop + 49
		elif svtype != 'TRA':
			pos_start = target_sv.start
			pos_stop = target_sv.stop
		else:
			pos_start = target_sv.start
			pos_stop = target_sv.start+1

		this_line = chrom+'\t'+str(pos_start)+'\t'+str(pos_stop)
		command1 = ('echo -e '+this_line).split(' ')
		command_tar = 'bedtools intersect -a -'.split(' ') + ('-b '+STR_file_TRF_target+' -f 0.5 -wa -wb').split(' ') 
		command_all = 'bedtools intersect -a -'.split(' ') + ('-b '+STR_file_TRF_all+' -f 0.5 -wa -wb').split(' ') 
		command_rms = 'bedtools intersect -a -'.split(' ') + ('-b '+STR_file_RM_Simpe+' -f 0.5 -wa -wb').split(' ') 
		command_mg = 'bedtools intersect -a -'.split(' ') + ('-b '+STR_file_MG+' -f 0.5 -wa -wb').split(' ') 
		#print('command1:', command1)
		#print('command_tar:', command_tar)
		#print('command_all:', command_all)
		#print('command_rms:', command_rms)
		#print('command_mg:', command_mg)

		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		str_tar_isecs = subprocess.run(command_tar, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		str_all_isecs = subprocess.run(command_all, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		str_rms_isecs = subprocess.run(command_rms, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		str_mg_isecs = subprocess.run(command_mg, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		#print('str_tar_isecs:', str_tar_isecs)
		#print('str_all_isecs:', str_all_isecs)
		#print('str_rms_isecs:', str_rms_isecs)
		#print('str_mg_isecs:', str_mg_isecs)

		str_tar_isecs = str_tar_isecs.split('\n')[:-1] # the last one is always an empty string
		str_all_isecs = str_all_isecs.split('\n')[:-1] # the last one is always an empty string
		str_rms_isecs = str_rms_isecs.split('\n')[:-1] # the last one is always an empty string
		str_mg_isecs = str_mg_isecs.split('\n')[:-1] # the last one is always an empty string
		#print('str_tar_isecs:', str_tar_isecs)
		#print('str_all_isecs:', str_all_isecs)
		#print('str_rms_isecs:', str_rms_isecs)
		#print('str_mg_isecs:', str_mg_isecs)

		if len(str_tar_isecs)>0:
			STR_bool = True
			str_start_list = []
			str_end_list = []
			period_len_list = []
			CN_list = []
			period_seq_list = []
			for str_isec in str_tar_isecs:
				_, _, _, str_isec_chrom, str_isec_start, str_isec_end, period_len, CN, period_seq = str_isec.split('\t')
				str_start_list.append(int(str_isec_start))
				str_end_list.append(int(str_isec_end))
				period_len_list.append(int(period_len))
				CN_list.append(float(CN))
				period_seq_list.append(period_seq)
			#print('str_start_list:', str_start_list)
			#print('str_end_list:', str_end_list)
			#print('period_len_list:', period_len_list)
			#print('CN_list:', CN_list)
			#print('period_seq_list:', period_seq_list)
			
		else:
			STR_bool = False

		rec.info['STR'] = STR_bool
		if STR_bool:
			rec.info['STR_REPEAT_LEN'] = ','.join([str(x) for x in period_len_list])
			rec.info['STR_REPEAT_SEQ'] = ','.join([str(x) for x in period_seq_list])
		if not STR_bool:
			if len(str_all_isecs)>0:
				rec.info['STR_TRF_OTHER'] = True
			if len(str_rms_isecs)>0:
				rec.info['STR_RM_SR'] = True
			if len(str_mg_isecs)>0:
				rec.info['STR_MG'] = True

		sample_supp_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()}
		visited_read_set = set()
		str_supp_dict_list = [{} for i in range(len(str_tar_isecs))]
		str_GT_list = []
		str_GQ_list = []

		for i_read, read in enumerate(fh_bam.fetch(chrom, pos_start-region_buffer_length, pos_stop+region_buffer_length)):
			if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
				if STR_bool and (svtype=='INS' or svtype=='DEL' or svtype=='DUP'):
					locus_read_name_list, num_repeat_list = str_signature(read, target_sv, str_start_list, str_end_list, period_len_list, CN_list, period_seq_list, k_s_dict, fa_handle, visited_read_set)
					visited_read_set.update([read.query_name])
					for i_str in range(len(str_tar_isecs)):
						if num_repeat_list[i_str] >= 0:
							read_name = locus_read_name_list[i_str]
							if read_name in str_supp_dict_list[i_str]:
								str_supp_dict_list[i_str][read_name] += num_repeat_list[i_str]
							else:
								str_supp_dict_list[i_str][read_name] = num_repeat_list[i_str]
				locus_read, CG_supp, SA_supp = sv_signiture(read, target_sv)
				sample_supp_dict['locus_reads'].update([locus_read])
				sample_supp_dict['CG_supp'].update([CG_supp])
				sample_supp_dict['SA_supp'].update([SA_supp])
		if STR_bool and (svtype=='INS' or svtype=='DEL' or svtype=='DUP'):
			for str_supp_dict in str_supp_dict_list:
				count_list = list(str_supp_dict.values())
				str_GT, str_GQ = infer_gt_str(count_list, p_err=0.05, svtype=svtype)
				str_GT_list.append(str_GT)
				str_GQ_list.append(str_GQ)
			rec.samples[0]['GT_STR'] = ','.join(str_GT_list)
			rec.samples[0]['GQ_STR'] = ','.join([str(x) for x in str_GQ_list])
			#print('str_num_repeats:', str_num_repeats)
			#print('[list(x.values()) for x in str_supp_dict_list]:')
			#print([list(x.values()) for x in str_supp_dict_list])
			#print('str_GT_list:', str_GT_list)
			#print('str_GQ_list:', str_GQ_list)
		sample_supp_dict['locus_reads'] -= {''}
		sample_supp_dict['CG_supp'] -= {''}
		sample_supp_dict['SA_supp'] -= {''}
		DV_s = len(sample_supp_dict['CG_supp'] | sample_supp_dict['SA_supp'])
		DR_s = len(sample_supp_dict['locus_reads']) - DV_s
		assert DR_s >= 0, 'problem with DR/DV, DR: '+str(DR_s)+', DV: '+str(DV_s)+', sv_id: '+str(sv_id)
		GT, GQ, p_11, p_01, p_00 = infer_gt(DR_s, DV_s, p_err)
		rec.samples[0]['RV'] = DV_s
		rec.samples[0]['RR'] = DR_s
		rec.samples[0]['GT_SV'] = GT
		rec.samples[0]['GQ_SV'] = GQ
		rec.samples[0]['P_11'] = p_11
		rec.samples[0]['P_01'] = p_01
		rec.samples[0]['P_00'] = p_00
		#print('DV_s:', DV_s, 'DR_s:', DR_s, 'GT:', GT, 'GQ:', GQ, 'p_00:', p_00, 'p_01:', p_01, 'p_11:', p_11)

		fh_vcf_out.write(rec)

	if verbose == 1:
		print('Finished all variants, count_skip_region:', count_skip_region)
		print('count_skip_sec:', count_skip_sec)

	fh_vcf_in.close()
	fh_vcf_out.close()
	fh_bam.close()


if __name__ == '__main__':

	arglist = sys.argv
	print('arglist:')
	print(arglist)

	vcf_in = arglist[1]
	vcf_out = arglist[2]
	bam_file = arglist[3]
	chrom = arglist[4]
	if chrom=='None':
		chrom = None
	n_sec = int(arglist[5])
	i_sec = int(arglist[6])

	make_VCF_GT(vcf_in, vcf_out, contig=chrom, bam_file=bam_file, n_sec=n_sec, i_sec=i_sec, verbose=1)
