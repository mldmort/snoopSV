import pysam
import glob
import pandas as pd
import numpy as np
import sys
from math import log10, factorial
import subprocess
from genSV import sv_class, infer_gt_sv, infer_gt_tr, sv_signiture, tr_signature

def make_VCF_GT(vcf_in, vcf_out, contig, bam_file, n_sec, i_sec, verbose=1):

	### genotyping setting
	mapping_quality_thr = 20
	region_buffer_length = 1000
	p_err = 0.05

	TR_file_TRF_target = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot_per-2-250.bed'
	TR_file_TRF_all = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed'
	TR_file_RM_Simpe = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Repeats_Masker_Simple_repeat.bed'
	TR_file_MG = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/hg38_ver13.bed'

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

	new_header_INFO = ['##INFO=<ID=TR,Number=0,Type=Flag,Description="TR region, target for TR genotyping">',
	'##INFO=<ID=TR_TRF_OTHER,Number=0,Type=Flag,Description="TR region, in TRF track but not a target for genotyping">',
	'##INFO=<ID=TR_RM_SR,Number=0,Type=Flag,Description="TR region, in repeat masker track, not a target for genotyping">',
	'##INFO=<ID=TR_MG,Number=0,Type=Flag,Description="TR region, in Gymreklab targets, but not in TR targets for genotyping">',
	'##INFO=<ID=TR_REPEAT_LEN,Number=.,Type=String,Description="TR repeat length">',
	'##INFO=<ID=TR_REPEAT_SEQ,Number=.,Type=String,Description="TR repeat sequence">'
	]
	new_header_FORMAT = ['##FORMAT=<ID=RV,Number=1,Type=Integer,Description="Number of reads supporting the variant sequence, from genSV">',
	'##FORMAT=<ID=RR,Number=1,Type=Integer,Description="Number of reads around the breakpoints supporting the reference sequence, from genSV">',
	'##FORMAT=<ID=GT_SV,Number=1,Type=String,Description="Genotype of the variant, from genSV">',
	'##FORMAT=<ID=GQ_SV,Number=1,Type=Integer,Description="Phred-scale genotype quality score, from genSV">',
	'##FORMAT=<ID=P_11,Number=1,Type=Float,Description="probability of 1/1 genotype, from genSV">',
	'##FORMAT=<ID=P_01,Number=1,Type=Float,Description="probability of 0/1 genotype, from genSV">',
	'##FORMAT=<ID=P_00,Number=1,Type=Float,Description="probability of 0/0 genotype, from genSV">',
	'##FORMAT=<ID=GT_TR_AL,Number=.,Type=String,Description="Genotype in TR region, from genSV">',
	'##FORMAT=<ID=GQ_TR_AL,Number=.,Type=String,Description="Phred-scale genotype quality score in TR region, from genSV">',
	'##FORMAT=<ID=GT_TR_LN,Number=.,Type=String,Description="Genotype in TR region, from genSV">',
	'##FORMAT=<ID=GQ_TR_LN,Number=.,Type=String,Description="Phred-scale genotype quality score in TR region, from genSV">'
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
		command_tar = 'bedtools intersect -a -'.split(' ') + ('-b '+TR_file_TRF_target+' -f 0.5 -wa -wb').split(' ') 
		command_all = 'bedtools intersect -a -'.split(' ') + ('-b '+TR_file_TRF_all+' -f 0.5 -wa -wb').split(' ') 
		command_rms = 'bedtools intersect -a -'.split(' ') + ('-b '+TR_file_RM_Simpe+' -f 0.5 -wa -wb').split(' ') 
		command_mg = 'bedtools intersect -a -'.split(' ') + ('-b '+TR_file_MG+' -f 0.5 -wa -wb').split(' ') 
		#print('command1:', command1)
		#print('command_tar:', command_tar)
		#print('command_all:', command_all)
		#print('command_rms:', command_rms)
		#print('command_mg:', command_mg)

		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		tr_tar_isecs = subprocess.run(command_tar, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		tr_all_isecs = subprocess.run(command_all, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		tr_rms_isecs = subprocess.run(command_rms, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		ps1 = subprocess.Popen(command1, stdout=subprocess.PIPE)
		tr_mg_isecs = subprocess.run(command_mg, stdin=ps1.stdout, check=True, capture_output=True, text=True).stdout
		#print('tr_tar_isecs:', tr_tar_isecs)
		#print('tr_all_isecs:', tr_all_isecs)
		#print('tr_rms_isecs:', tr_rms_isecs)
		#print('tr_mg_isecs:', tr_mg_isecs)

		tr_tar_isecs = tr_tar_isecs.split('\n')[:-1] # the last one is always an empty string
		tr_all_isecs = tr_all_isecs.split('\n')[:-1] # the last one is always an empty string
		tr_rms_isecs = tr_rms_isecs.split('\n')[:-1] # the last one is always an empty string
		tr_mg_isecs = tr_mg_isecs.split('\n')[:-1] # the last one is always an empty string
		#print('tr_tar_isecs:', tr_tar_isecs)
		#print('tr_all_isecs:', tr_all_isecs)
		#print('tr_rms_isecs:', tr_rms_isecs)
		#print('tr_mg_isecs:', tr_mg_isecs)

		if len(tr_tar_isecs)>0:
			TR_bool = True
			tr_start_list = []
			tr_end_list = []
			period_len_list = []
			CN_list = []
			period_seq_list = []
			for tr_isec in tr_tar_isecs:
				_, _, _, tr_isec_chrom, tr_isec_start, tr_isec_end, period_len, CN, period_seq = tr_isec.split('\t')
				tr_start_list.append(int(tr_isec_start))
				tr_end_list.append(int(tr_isec_end))
				period_len_list.append(int(period_len))
				CN_list.append(float(CN))
				period_seq_list.append(period_seq)
			#print('tr_start_list:', tr_start_list)
			#print('tr_end_list:', tr_end_list)
			#print('period_len_list:', period_len_list)
			#print('CN_list:', CN_list)
			#print('period_seq_list:', period_seq_list)
		else:
			TR_bool = False

		rec.info['TR'] = TR_bool
		if TR_bool:
			rec.info['TR_REPEAT_LEN'] = ','.join([str(x) for x in period_len_list])
			rec.info['TR_REPEAT_SEQ'] = ','.join([str(x) for x in period_seq_list])
		if not TR_bool:
			if len(tr_all_isecs)>0:
				rec.info['TR_TRF_OTHER'] = True
			if len(tr_rms_isecs)>0:
				rec.info['TR_RM_SR'] = True
			if len(tr_mg_isecs)>0:
				rec.info['TR_MG'] = True

		sample_supp_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()}
		visited_read_set = set()
		tr_supp_al_list = [[] for i in range(len(tr_tar_isecs))]
		tr_supp_ln_list = [[] for i in range(len(tr_tar_isecs))]
		tr_GT_al_list = []
		tr_GT_ln_list = []
		tr_GQ_al_list = []
		tr_GQ_ln_list = []

		for i_read, read in enumerate(fh_bam.fetch(chrom, pos_start-region_buffer_length, pos_stop+region_buffer_length)):
			if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
				if TR_bool and (svtype=='INS' or svtype=='DEL'):
					locus_read_name_list, num_repeat_al_list, num_repeat_ln_list = tr_signature(read, target_sv, tr_start_list, tr_end_list, period_len_list, CN_list, period_seq_list, k_s_dict, fa_handle, visited_read_set)
					visited_read_set.update([read.query_name])
					for i_tr in range(len(tr_tar_isecs)):
						if num_repeat_al_list[i_tr] >= 0:
							tr_supp_al_list[i_tr].append(num_repeat_al_list[i_tr])
						if num_repeat_ln_list[i_tr] >= 0:
							tr_supp_ln_list[i_tr].append(num_repeat_ln_list[i_tr])
				locus_read, CG_supp, SA_supp = sv_signiture(read, target_sv)
				sample_supp_dict['locus_reads'].update([locus_read])
				sample_supp_dict['CG_supp'].update([CG_supp])
				sample_supp_dict['SA_supp'].update([SA_supp])
		if TR_bool and (svtype=='INS' or svtype=='DEL'):
			for count_list in tr_supp_al_list:
				tr_GT, tr_GQ = infer_gt_tr(count_list, p_err=0.05, svtype=svtype)
				tr_GT_al_list.append(tr_GT)
				tr_GQ_al_list.append(tr_GQ)
			for count_list in tr_supp_ln_list:
				tr_GT, tr_GQ = infer_gt_tr(count_list, p_err=0.05, svtype=svtype)
				tr_GT_ln_list.append(tr_GT)
				tr_GQ_ln_list.append(tr_GQ)
			rec.samples[0]['GT_TR_AL'] = ','.join(tr_GT_al_list)
			rec.samples[0]['GQ_TR_AL'] = ','.join([str(x) for x in tr_GQ_al_list])
			rec.samples[0]['GT_TR_LN'] = ','.join(tr_GT_ln_list)
			rec.samples[0]['GQ_TR_LN'] = ','.join([str(x) for x in tr_GQ_ln_list])
		sample_supp_dict['locus_reads'] -= {''}
		sample_supp_dict['CG_supp'] -= {''}
		sample_supp_dict['SA_supp'] -= {''}
		DV_s = len(sample_supp_dict['CG_supp'] | sample_supp_dict['SA_supp'])
		DR_s = len(sample_supp_dict['locus_reads']) - DV_s
		assert DR_s >= 0, 'problem with DR/DV, DR: '+str(DR_s)+', DV: '+str(DV_s)+', sv_id: '+str(sv_id)
		GT, GQ, p_11, p_01, p_00 = infer_gt_sv(DR_s, DV_s, p_err)
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
