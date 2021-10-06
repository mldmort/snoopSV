import pysam
import glob
import pandas as pd
import numpy as np
import sys
from math import log10, factorial
import subprocess
from genSV import sv_class, infer_gt_sv, infer_gt_tr, infer_gt_tr_phased, sv_signiture, tr_signature, tr_signature_2
import argparse
import pickle

### global variables
new_header_INFO = [
'##INFO=<ID=SKIP_REGION,Number=0,Type=Flag,Description="Skip this call for genotyping, since it falls into skipped regions">',
'##INFO=<ID=SKIP_TR,Number=0,Type=Flag,Description="Skip this call for genotyping, since it falls into TR regions">',
'##INFO=<ID=TR_ANNOT,Number=0,Type=Flag,Description="Tandem Repeat annotation genotyped, from genSV">',
'##INFO=<ID=TR,Number=0,Type=Flag,Description="TR region, target for TR genotyping">',
'##INFO=<ID=TR_TRF_OTHER,Number=0,Type=Flag,Description="TR region, in TRF track but not a target for genotyping">',
'##INFO=<ID=TR_RM_SR,Number=0,Type=Flag,Description="TR region, in repeat masker track, not a target for genotyping">',
'##INFO=<ID=TR_MG,Number=0,Type=Flag,Description="TR region, in Gymreklab targets, but not in TR targets for genotyping">',
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
'##FORMAT=<ID=CN_TR_BP_H1,Number=.,Type=String,Description="Counts of base pairs in TR region for H1, from genSV">',
'##FORMAT=<ID=CN_TR_BP_H2,Number=.,Type=String,Description="Counts of base pairs in TR region for H2, from genSV">',
'##FORMAT=<ID=CN_TR_BP_H0,Number=.,Type=String,Description="Counts of base pairs in TR region for H0, from genSV">'
]
#TR_file_TRF_target = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot_per-2-250.bed'
TR_file_TRF_target = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed'
TR_file_TRF_all = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Simple_Repeats_TRF_annot.bed'
TR_file_RM_Simpe = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/Repeats_Masker_Simple_repeat.bed'
TR_file_MG = '/home/smmortazavi/HUMAN_DATA/REF/REPEATS/hg38_ver13.bed'

skip_region_list = [\
	{'chrom':'chr1', 'start':143150000, 'stop':149900000}, \
	{'chrom':'chr16', 'start':46380000, 'stop':46425000}, \
	{'chrom':'chr1', 'start':125060000, 'stop':125200000} \
]

def skip_region(skip_region_list, chrom, start, stop):

	if (chrom==skip_region_list[0]['chrom']) and ((start > skip_region_list[0]['start']) and (start < skip_region_list[0]['stop'])) and ((stop > skip_region_list[0]['start']) and (stop < skip_region_list[0]['stop'])):
		return True

	if (chrom==skip_region_list[1]['chrom']) and ((start > skip_region_list[1]['start']) and (start < skip_region_list[1]['stop'])) and ((stop > skip_region_list[1]['start']) and (stop < skip_region_list[1]['stop'])):
		return True

	if (chrom==skip_region_list[2]['chrom']) and ((start > skip_region_list[2]['start']) and (start < skip_region_list[2]['stop'])) and ((stop > skip_region_list[2]['start']) and (stop < skip_region_list[2]['stop'])):
		return True

	return False

def get_phased_gt(GT_SV, RV_P, RV_M):
	if GT_SV == './.':
		return GT_SV
	assert ((RV_P>=0) and (RV_M>=0)), 'problem RV_P: '+str(RV_P)+' ,RV_M: '+str(RV_M)
	assert (GT_SV=='1/1' or GT_SV=='0/1' or GT_SV=='0/0'), 'problem GT_SV: '+GT_SV
	if GT_SV == '1/1':
		if (RV_P > 0) and (RV_M > 0):
			GT_SV_PH = '1|1'
		elif (RV_P == 0) or (RV_M == 0):
			GT_SV_PH = GT_SV
		else:
			assert 0==1, 'problem in get_phased_gt, GT_SV: '+GT_SV+' ,RV_P: '+str(RV_P)+' ,RV_M: '+str(RV_M)
	elif GT_SV == '0/1':
		if (RV_P > 0) and (RV_M > 0):
			GT_SV_PH = '1|1'
		elif (RV_P > 0) and (RV_M == 0):
			GT_SV_PH = '1|0'
		elif (RV_P == 0) and (RV_M > 0):
			GT_SV_PH = '0|1'
		elif (RV_P == 0) and (RV_M == 0):
			GT_SV_PH = GT_SV
		else:
			assert 0==1, 'problem in get_phased_gt, GT_SV: '+GT_SV+' ,RV_P: '+str(RV_P)+' ,RV_M: '+str(RV_M)
	elif GT_SV == '0/0':
		if (RV_P > 0) and (RV_M > 0):
			GT_SV_PH = '1|1'
		elif (RV_P > 0) and (RV_M == 0):
			GT_SV_PH = '1|0'
		elif (RV_P == 0) and (RV_M > 0):
			GT_SV_PH = '0|1'
		elif (RV_P == 0) and (RV_M == 0):
			GT_SV_PH = GT_SV
		else:
			assert 0==1, 'problem in get_phased_gt, GT_SV: '+GT_SV+' ,RV_P: '+str(RV_P)+' ,RV_M: '+str(RV_M)

	return GT_SV_PH

def GT_TR(tr_annot_file, vcf_in, vcf_out, contig, sample_bam_file, n_sec, i_sec, tr_span_max, verbose=1):

	### genotyping setting
	mapping_quality_thr = 20
	region_buffer_length = 1000

	sample_bam_dict = {}
	with open(sample_bam_file, 'r') as fh:
		for line in fh.readlines():
			line = line.strip().split()
			sample_bam_dict[line[0]] = line[1]

	k_s_dict = {}
	#**p_val_thr = 0.07
	#**for k in range(2,262):
	#**	for s_sel in range(k,0,-1):
	#**		p_val = factorial(k)/factorial(s_sel)/factorial(k-s_sel)/float(4**s_sel)
	#**		if p_val > p_val_thr:
	#**			k_s_dict[k] = s_sel + 1
	#**			break

	if contig:
		command = ['awk', 'BEGIN{FS="\t";OFS="\t"}$1=="'+contig+'"', tr_annot_file]
	else:
		command = ['awk', 'BEGIN{FS="\t";OFS="\t"}$1!~"#"', tr_annot_file]
	ps = subprocess.Popen(command, stdout=subprocess.PIPE)
	n_calls = int(subprocess.run(['wc', '-l'], stdin=ps.stdout, check=True, capture_output=True, text=True).stdout)

	### i_sec should go from 0 to n_sec to cover all calls
	i_rec_start = i_sec * int(n_calls / n_sec)
	i_rec_end = (i_sec + 1) * int(n_calls / n_sec)
	if (i_sec == n_sec):
		i_rec_end = n_calls
	if verbose == 1:
		print('n_calls:',  n_calls)
		print('i_sec:', i_sec)
		print('n_sec:', n_sec)
		print('i_rec_start:', i_rec_start)
		print('i_rec_end:', i_rec_end)
	ps1 = subprocess.Popen(command, stdout=subprocess.PIPE)
	ps2 = subprocess.Popen(['head', '-n', str(i_rec_end)], stdin=ps1.stdout, stdout=subprocess.PIPE)
	out = subprocess.run(['tail', '-n', str(i_rec_end-i_rec_start)], stdin=ps2.stdout, check=True, capture_output=True, text=True).stdout
	tr_annots_list = out.split('\n')[:-1] # the last one is always an empty string
	#print('tr_annots_list[0]:', tr_annots_list[0])
	#print('tr_annots_list[-1]:', tr_annots_list[-1])

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	for line in new_header_INFO+new_header_FORMAT:
		header_in.add_line(line)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	count_skip_region = 0
	for tr_annot in tr_annots_list:
		tr_chrom, tr_start, tr_end, period_len, CN, period_seq = tr_annot.split('\t')
		tr_start = int(tr_start)
		tr_end = int(tr_end)
		period_len = int(period_len)
		CN = float(CN)

		rec = fh_vcf_out.new_record(contig=tr_chrom, start=tr_start-1, stop=tr_end, alleles=('.', '.'))

		rec.info['TR_REPEAT_LEN'] = str(period_len)
		rec.info['TR_REPEAT_SEQ'] = period_seq
		rec.info['TR_REPEAT_START'] = str(tr_start)
		rec.info['TR_REPEAT_END'] = str(tr_end)
		rec.info['TR_REPEAT_CN'] = str(CN)
		rec.info['TR_ANNOT'] = True

		if skip_region(skip_region_list, tr_chrom, tr_start, tr_end):
			count_skip_region += 1
			rec.info['SKIP_REGION'] = True
			fh_vcf_out.write(rec)
			continue

		for sample, bam_file in sample_bam_dict.items():
			fh_bam = pysam.AlignmentFile(bam_file, 'rb')
			visited_read_set = set()
			tr_supp_al = {'H1':[], 'H2':[], 'H0':[]}
			tr_supp_ln = {'H1':[], 'H2':[], 'H0':[]}
			tr_supp_bp = {'H1':[], 'H2':[], 'H0':[]}

			for i_read, read in enumerate(fh_bam.fetch(tr_chrom, max(0,tr_start-region_buffer_length), tr_end+region_buffer_length)):
				if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
					locus_read_name, num_repeat_al, num_repeat_ln, num_bp = tr_signature_2(read, tr_start, tr_end, period_len, CN, period_seq, k_s_dict, visited_read_set)
					visited_read_set.update([locus_read_name])
					if num_repeat_al >= 0:
						if read.has_tag('HP'):
							HP = read.get_tag(tag='HP')
							tr_supp_al['H'+str(HP)].append(num_repeat_al)
						else:
							tr_supp_al['H0'].append(num_repeat_al)
					if num_repeat_ln >= 0:
						if read.has_tag('HP'):
							HP = read.get_tag(tag='HP')
							tr_supp_bp['H'+str(HP)].append(num_bp)
							tr_supp_ln['H'+str(HP)].append(num_repeat_ln)
						else:
							tr_supp_bp['H0'].append(num_bp)
							tr_supp_ln['H0'].append(num_repeat_ln)
			fh_bam.close()

			tr_GT_ln = infer_gt_tr_phased(tr_supp_ln)
			tr_GT_bp = infer_gt_tr_phased(tr_supp_bp)
			rec.samples[sample]['GT_TR_LN'] = tr_GT_ln
			rec.samples[sample]['GT_TR_BP'] = tr_GT_bp

			temp = '|'.join([str(x) for x in tr_supp_ln['H1']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_LN_H1'] = temp

			temp = '|'.join([str(x) for x in tr_supp_ln['H2']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_LN_H2'] = temp

			temp = '|'.join([str(x) for x in tr_supp_ln['H0']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_LN_H0'] = temp

			temp = '|'.join([str(x) for x in tr_supp_bp['H1']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_BP_H1'] = temp

			temp = '|'.join([str(x) for x in tr_supp_bp['H2']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_BP_H2'] = temp

			temp = '|'.join([str(x) for x in tr_supp_bp['H0']])
			if temp == '':
				temp = '.'
			rec.samples[sample]['CN_TR_BP_H0'] = temp

		fh_vcf_out.write(rec)

	if verbose == 1:
		print('Finished all variants')
		print('count_skip_region:', count_skip_region)

	fh_vcf_out.close()

def GT_nonTR(vcf_in, vcf_out, contig, sample_bam_file, n_sec, i_sec, verbose=1):

	### genotyping setting
	mapping_quality_thr = 20
	region_buffer_length = 1000
	SV_p_err = 0.01

	sample_bam_dict = {}
	with open(sample_bam_file, 'r') as fh:
		for line in fh.readlines():
			line = line.strip().split()
			sample_bam_dict[line[0]] = line[1]

	if contig:
		command = ('bcftools query -r '+contig+' -f %CHROM\\n '+vcf_in).split()
	else:
		command = ('bcftools query -f %CHROM\\n '+vcf_in).split()
	ps = subprocess.Popen(command, stdout=subprocess.PIPE)
	n_calls = int(subprocess.run(['wc', '-l'], stdin=ps.stdout, check=True, capture_output=True, text=True).stdout)
	### i_sec should go from 0 to n_sec to cover all calls
	i_rec_start = i_sec * int(n_calls / n_sec)
	i_rec_end = (i_sec + 1) * int(n_calls / n_sec)
	if (i_sec == n_sec):
		i_rec_end = n_calls
	if verbose == 1:
		print('n_calls:',  n_calls)
		print('i_sec:', i_sec)
		print('n_sec:', n_sec)
		print('i_rec_start:', i_rec_start)
		print('i_rec_end:', i_rec_end)

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	for line in new_header_INFO+new_header_FORMAT:
		header_in.add_line(line)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	count_skip_region = 0
	count_skip_sec = 0
	count_skip_tr = 0
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

		if skip_region(skip_region_list, chrom, start, stop):
			count_skip_region += 1
			rec.info['SKIP_REGION'] = True
			fh_vcf_out.write(rec)
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
			rec.info['TR_REPEAT_START'] = ','.join([str(x) for x in tr_start_list])
			rec.info['TR_REPEAT_END'] = ','.join([str(x) for x in tr_end_list])
			rec.info['TR_REPEAT_CN'] = ','.join([str(x) for x in CN_list])
		if not TR_bool:
			if len(tr_all_isecs)>0:
				rec.info['TR_TRF_OTHER'] = True
			if len(tr_rms_isecs)>0:
				rec.info['TR_RM_SR'] = True
			if len(tr_mg_isecs)>0:
				rec.info['TR_MG'] = True

		if TR_bool and (svtype=='INS' or svtype=='DEL'):
			count_skip_tr += 1
			rec.info['SKIP_TR'] = True
			fh_vcf_out.write(rec)
			continue

		for sample, bam_file in sample_bam_dict.items():
			fh_bam = pysam.AlignmentFile(bam_file, 'rb')
			read_supp_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()}
			read_supp_P_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### paternal
			read_supp_M_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### maternal
			read_supp_N_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### not phased

			for i_read, read in enumerate(fh_bam.fetch(chrom, max(0,pos_start-region_buffer_length), pos_stop+region_buffer_length)):
				if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
					locus_read, CG_supp, SA_supp = sv_signiture(read, target_sv)
					read_supp_dict['locus_reads'].update([locus_read])
					read_supp_dict['CG_supp'].update([CG_supp])
					read_supp_dict['SA_supp'].update([SA_supp])
					if read.has_tag('HP'):
						HP = read.get_tag(tag='HP')
						if HP == 1:
							read_supp_P_dict['locus_reads'].update([locus_read])
							read_supp_P_dict['CG_supp'].update([CG_supp])
							read_supp_P_dict['SA_supp'].update([SA_supp])
						elif HP == 2:
							read_supp_M_dict['locus_reads'].update([locus_read])
							read_supp_M_dict['CG_supp'].update([CG_supp])
							read_supp_M_dict['SA_supp'].update([SA_supp])
						else:
							assert 0==1, 'problem with HP in read: ' + read.query_name
					else:
						read_supp_N_dict['locus_reads'].update([locus_read])
						read_supp_N_dict['CG_supp'].update([CG_supp])
						read_supp_N_dict['SA_supp'].update([SA_supp])
			fh_bam.close()
			read_supp_dict['locus_reads'] -= {''}
			read_supp_dict['CG_supp'] -= {''}
			read_supp_dict['SA_supp'] -= {''}
			read_supp_P_dict['locus_reads'] -= {''}
			read_supp_P_dict['CG_supp'] -= {''}
			read_supp_P_dict['SA_supp'] -= {''}
			read_supp_M_dict['locus_reads'] -= {''}
			read_supp_M_dict['CG_supp'] -= {''}
			read_supp_M_dict['SA_supp'] -= {''}
			read_supp_N_dict['locus_reads'] -= {''}
			read_supp_N_dict['CG_supp'] -= {''}
			read_supp_N_dict['SA_supp'] -= {''}
			DV_s = len(read_supp_dict['CG_supp'] | read_supp_dict['SA_supp'])
			DR_s = len(read_supp_dict['locus_reads']) - DV_s
			assert DR_s >= 0, 'problem with DR/DV, DR: '+str(DR_s)+', DV: '+str(DV_s)+', sv_id: '+str(sv_id)
			DV_s_P = len(read_supp_P_dict['CG_supp'] | read_supp_P_dict['SA_supp'])
			DR_s_P = len(read_supp_P_dict['locus_reads']) - DV_s_P
			assert DR_s_P >= 0, 'problem with P DR/DV, DR: '+str(DR_s_P)+', DV: '+str(DV_s_P)+', sv_id: '+str(sv_id)
			DV_s_M = len(read_supp_M_dict['CG_supp'] | read_supp_M_dict['SA_supp'])
			DR_s_M = len(read_supp_M_dict['locus_reads']) - DV_s_M
			assert DR_s_M >= 0, 'problem with M DR/DV, DR: '+str(DR_s_M)+', DV: '+str(DV_s_M)+', sv_id: '+str(sv_id)
			DV_s_N = len(read_supp_N_dict['CG_supp'] | read_supp_N_dict['SA_supp'])
			DR_s_N = len(read_supp_N_dict['locus_reads']) - DV_s_N
			assert DR_s_N >= 0, 'problem with N DR/DV, DR: '+str(DR_s_N)+', DV: '+str(DV_s_N)+', sv_id: '+str(sv_id)
			assert DV_s_P+DV_s_M+DV_s_N == DV_s, 'problem with P/M/N DV_s, DV_s_P: '+str(DV_s_P)+', DV_s_M: '+str(DV_s_M)+', DV_s_N: '+str(DV_s_N)+', DV_s: '+str(DV_s)
			assert DR_s_P+DR_s_M+DR_s_N == DR_s, 'problem with P/M/N DR_s, DR_s_P: '+str(DR_s_P)+', DR_s_M: '+str(DR_s_M)+', DR_s_N: '+str(DR_s_N)+', DR_s: '+str(DR_s)
			GT, GQ, p_11, p_01, p_00, SQ = infer_gt_sv(DR_s, DV_s, p_err=SV_p_err)
			GT_PH = get_phased_gt(GT, DV_s_P, DV_s_M)
			rec.samples[sample]['RV'] = DV_s
			rec.samples[sample]['RR'] = DR_s
			rec.samples[sample]['RV_P'] = DV_s_P
			rec.samples[sample]['RR_P'] = DR_s_P
			rec.samples[sample]['RV_M'] = DV_s_M
			rec.samples[sample]['RR_M'] = DR_s_M
			rec.samples[sample]['RV_N'] = DV_s_N
			rec.samples[sample]['RR_N'] = DR_s_N
			rec.samples[sample]['GT_SV'] = GT
			rec.samples[sample]['GQ_SV'] = GQ
			rec.samples[sample]['P_11'] = p_11
			rec.samples[sample]['P_01'] = p_01
			rec.samples[sample]['P_00'] = p_00
			rec.samples[sample]['SQ_SV'] = SQ
			rec.samples[sample]['GT_SV_PH'] = GT_PH
			#print('DV_s:', DV_s, 'DR_s:', DR_s, 'GT:', GT, 'GQ:', GQ, 'p_00:', p_00, 'p_01:', p_01, 'p_11:', p_11)

		fh_vcf_out.write(rec)

	if verbose == 1:
		print('Finished all variants')
		print('count_skip_region:', count_skip_region)
		print('count_skip_sec:', count_skip_sec)
		print('count_skip_tr:', count_skip_tr)

	fh_vcf_in.close()
	fh_vcf_out.close()

def get_data(files, columns=None, header=None, index_col=None):

	data = pd.read_table(files[0], header=header, index_col=index_col, sep='\t', keep_default_na=False)
	for f in files[1:]:
		data = pd.concat([data, pd.read_table(f, header=None, index_col=None, sep='\t', keep_default_na=False)], ignore_index=True, axis=0)
	if columns != None:
		data.columns = columns
	return data

def get_data_lin(annot_in, cov_in):

	sample_cov_dict = {}
	with open(cov_in, 'r') as fh:
		for line in fh.readlines():
			sample = line.strip().split('\t')[0]
			cov_file = line.strip().split('\t')[1]
			sample_cov_dict[sample] = get_data([cov_file], header=0, index_col=0)

	main_chroms_list = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']

	data_annot = pd.read_table(annot_in, sep='\t', header=0)
	#print(data_annot)
	data_annot = data_annot.loc[data_annot.chr2.isin(main_chroms_list)]
	#print(data_annot)

	if data_annot.empty:
		return pd.DataFrame()

	### data_lin will have columns_common columns plus these ones: rv, dp, dp_fr, mapq, sample
	columns_common = ['chrom', 'pos', 'end', 'id', 'svtype', 'chr2', 'end2']

	data_lin = pd.DataFrame()
	for index in data_annot.index.tolist():
		svtype = data_annot.loc[index, 'svtype']
		samples_list = str(data_annot.loc[index, 'samples']).split(',')
		mapq_list = str(data_annot.loc[index, 'mapq_all_samples']).split(',')
		dp_fr_list = str(data_annot.loc[index, 'dp_fr_samples']).split(',')
		dp_list = str(data_annot.loc[index, 'dp_samples']).split(',')
		rv_list = str(data_annot.loc[index, 'rv']).split(',')
		for i_sam, sample in enumerate(samples_list):
			mapq = mapq_list[i_sam]
			if svtype!='TRA':
				dp_fr = dp_fr_list[i_sam]
			else:
				dp_fr = '.'
			dp = dp_list[i_sam]
			rv = rv_list[i_sam]
			temp = pd.DataFrame(data_annot.loc[[index], columns_common])
			temp['rv'] = rv
			temp['dp'] = dp
			temp['dp_fr'] = dp_fr
			temp['mapq'] = mapq
			temp['sample'] = sample
			if svtype=='TRA':
				temp = pd.concat([temp, temp], axis=0, ignore_index=True)
				temp.loc[0, 'svtype'] = 'TRA_1'
				temp.loc[1, 'svtype'] = 'TRA_2'
				temp.loc[0, 'dp'] = temp.loc[0, 'dp'].split(':')[0]
				temp.loc[1, 'dp'] = temp.loc[1, 'dp'].split(':')[1]
			data_lin = pd.concat([data_lin, temp], axis=0, ignore_index=True)

	data_lin['rv'] = data_lin['rv'].astype(float)
	data_lin['dp'] = data_lin['dp'].astype(float)
	data_lin['mapq'] = data_lin['mapq'].astype(float)

	sub_data = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2')]
	if not sub_data.empty:
		### split dp_fr
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_1'] = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr'].str.split(pat=':', expand=True)[0]
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_2'] = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr'].str.split(pat=':', expand=True)[1]
		### correct depth
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_cor'] = (data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_1'].astype(float) + data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_2'].astype(float) + 1e-7) / 2.

	sub_data = data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2')]
	if not sub_data.empty:
		### correct depth
		data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2'), 'dp_cor'] = data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2'), 'dp'] + 1e-7

	### compute DV/depth with corrected depth
	data_lin['rv/depth'] = data_lin.rv / data_lin.dp_cor

	### compute depth/cov
	for key in sample_cov_dict:
		for chrom in ['chr'+str(n) for n in range(1,23)]+['chrX', 'chrY']:
			data_lin.loc[(data_lin['sample']==key) & (data_lin.chrom==chrom) & (data_lin.svtype!='TRA_2'), 'cov'] = sample_cov_dict[key].loc[chrom, 'mean']
			data_lin.loc[(data_lin['sample']==key) & (data_lin.chr2==chrom) & (data_lin.svtype=='TRA_2'), 'cov'] = sample_cov_dict[key].loc[chrom, 'mean']
	data_lin['depth/cov'] = data_lin['dp_cor'] / data_lin['cov']

	return data_lin

def apply_ML_models(data_lin, model_labels, model_files, model_vars):

	for i_model, model_file in enumerate(model_files):

		model_var_list = model_vars[i_model]
		model_label = model_labels[i_model]
		print('working on model:', model_label)

		with open(model_file, 'rb') as fh:
			clf = pickle.load(fh)

		X_test = np.array(data_lin[model_var_list])
		predict_proba = clf.predict_proba(X_test)[:,1] # [:,0]: probablity of class 0; [:,1]: probablity of class 1

		data_lin[model_label] = predict_proba

	return data_lin

def write_vcf(vcf_in, data_lin, model_labels, vcf_out):

	new_headers = []
	for model_label in model_labels:
		new_headers.append('##FORMAT=<ID='+model_label+',Number=1,Type=String,Description="ML model score">')

	new_headers.append('##FORMAT=<ID=DP,Number=1,Type=String,Description="regional read depth">')
	new_headers.append('##FORMAT=<ID=COV,Number=1,Type=String,Description="average chromosome coverage">')
	new_headers.append('##FORMAT=<ID=MapQ,Number=1,Type=String,Description="average mapping quality of all intersecting reads">')
	new_headers.append('##FORMAT=<ID=DP_o_COV,Number=1,Type=String,Description="regional read depth over chrom coverage">')
	new_headers.append('##FORMAT=<ID=RV_o_DP,Number=1,Type=String,Description="read support over regional read depth">')
	new_headers.append('##FORMAT=<ID=DP_FR,Number=1,Type=String,Description="read depth of flanking regions">')

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	for line in new_headers:
		header_in.add_line(line)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	for i_rec, rec in enumerate(fh_vcf_in.fetch()):
		sv_id = rec.id
		svtype = rec.info['SVTYPE']

		if not data_lin.empty:
			if svtype != 'TRA':
				sub_data = data_lin.loc[data_lin.id==sv_id]
				for index in sub_data.index.tolist():
					sample = sub_data.loc[index, 'sample']
					dp = sub_data.loc[index, 'dp_cor']
					dp_fr = '|'.join(sub_data.loc[index, 'dp_fr'].split(':'))
					cov = sub_data.loc[index, 'cov']
					dp_o_cov = sub_data.loc[index, 'depth/cov']
					mapq = sub_data.loc[index, 'mapq']
					rv_o_dp = sub_data.loc[index, 'rv/depth']
					model_scores = {}
					for model in model_labels:
						model_scores[model] = sub_data.loc[index, model]

					rec.samples[sample]['DP'] = '{:.2f}'.format(dp)
					rec.samples[sample]['COV'] = '{:.2f}'.format(cov)
					rec.samples[sample]['MapQ'] = '{:.2f}'.format(mapq)
					rec.samples[sample]['DP_o_COV'] = '{:.2f}'.format(dp_o_cov)
					rec.samples[sample]['RV_o_DP'] = '{:.2f}'.format(rv_o_dp)
					rec.samples[sample]['DP_FR'] = str(dp_fr)
					for model, score in model_scores.items():
						rec.samples[sample][model] = '{:.2f}'.format(score)
			else:
				sub_data_1 = data_lin.loc[(data_lin.id==sv_id) & (data_lin.svtype=='TRA_1')]
				sub_data_2 = data_lin.loc[(data_lin.id==sv_id) & (data_lin.svtype=='TRA_2')]
				#print('sub_data_1:')
				#print(sub_data_1)
				#print('sub_data_2:')
				#print(sub_data_2)
				for index in sub_data_1.index.tolist():
					sample = sub_data_1.loc[index, 'sample']
					mapq = sub_data_1.loc[index, 'mapq']
					dp_1 = sub_data_1.loc[index, 'dp_cor']
					dp_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'dp_cor'].values[0]
					cov_1 = sub_data_1.loc[index, 'cov']
					cov_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'cov'].values[0]
					dp_o_cov_1 = sub_data_1.loc[index, 'depth/cov']
					dp_o_cov_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'depth/cov'].values[0]
					rv_o_dp_1 = sub_data_1.loc[index, 'rv/depth']
					rv_o_dp_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'rv/depth'].values[0]
					model_scores_1 = {}
					model_scores_2 = {}
					for model in model_labels:
						model_scores_1[model] = sub_data_1.loc[index, model]
						model_scores_2[model] = sub_data_2.loc[sub_data_2['sample']==sample, model].values[0]

					rec.samples[sample]['DP'] = '{:.2f}'.format(dp_1)+'|'+'{:.2f}'.format(dp_2)
					rec.samples[sample]['COV'] = '{:.2f}'.format(cov_1)+'|'+'{:.2f}'.format(cov_2)
					rec.samples[sample]['MapQ'] = '{:.2f}'.format(mapq)
					rec.samples[sample]['DP_o_COV'] = '{:.2f}'.format(dp_o_cov_1)+'|'+'{:.2f}'.format(dp_o_cov_2)
					rec.samples[sample]['RV_o_DP'] = '{:.2f}'.format(rv_o_dp_1)+'|'+'{:.2f}'.format(rv_o_dp_2)
					for model, score_1 in model_scores_1.items():
						score_2 = model_scores_2[model]
						rec.samples[sample][model] = '{:.2f}'.format(score_1)+'|'+'{:.2f}'.format(score_2)
			
		fh_vcf_out.write(rec)

	fh_vcf_in.close()
	fh_vcf_out.close()
			
	
def SCORE_VCF(vcf_in, annot_in, cov_in, models, vcf_out):

	data_lin = get_data_lin(annot_in, cov_in)
	#print('data_lin:')
	#print(data_lin)

	model_labels = []
	model_files = []
	model_vars = []
	with open(models, 'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')
			model_labels.append(line[0])
			model_files.append(line[1])
			model_vars.append(line[2].split(','))
	#print('model_labels:', model_labels)
	#print('model_files:', model_files)
	#print('model_vars:', model_vars)
	if not data_lin.empty:
		data_lin = apply_ML_models(data_lin, model_labels, model_files, model_vars)
	#print('data_lin:')
	#print(data_lin)

	write_vcf(vcf_in, data_lin, model_labels, vcf_out)


def run_nontr(args):
	vcf_in = args.vcf_in
	vcf_out = args.vcf_out
	sample_bam_file = args.sample_bam_file
	chrom = args.contig
	n_sec = args.n_section
	i_sec = args.i_section
	for x, y in args.__dict__.items():
		print(x,':', y)

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
	for x, y in args.__dict__.items():
		print(x,':', y)

	GT_TR(tr_annot_file=tr_annot_file, vcf_in=vcf_in, vcf_out=vcf_out, contig=chrom, sample_bam_file=sample_bam_file, n_sec=n_sec, i_sec=i_sec, tr_span_max=tr_span_max, verbose=1)

def run_score(args):
	vcf_in = args.vcf_in
	annot_in = args.annot_in
	cov_in = args.coverage_file
	vcf_out = args.vcf_out
	models = args.models
	for x, y in args.__dict__.items():
		print(x,':', y)

	SCORE_VCF(vcf_in=vcf_in, annot_in=annot_in, cov_in=cov_in, models=models, vcf_out=vcf_out)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description='genSV is a genotyper for cohort SV calls from third generation data (ONT/PacBio)')
	subparsers = parser.add_subparsers(help='available sub-commands')

	parser_nontr = subparsers.add_parser('nontr', help='process non-TR SVs')
	parser_nontr.add_argument('-v', '--vcf_in', required=True, help='input VCF file')
	parser_nontr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_nontr.add_argument('-s', '--sample_bam_file', required=True, help='a map between samples and bam files as a tab delimited text file. First column is the samples, and second column is the absolute path to the bam files')
	parser_nontr.add_argument('-c', '--contig', default=None, help='contig name. If used the input VCF should be indexed')
	parser_nontr.add_argument('-n', '--n_section', type=int, default=1, help='number of sections in the input VCF file')
	parser_nontr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
	parser_nontr.set_defaults(func=run_nontr)

	parser_tr = subparsers.add_parser('tr', help='process TR annotations')
	parser_tr.add_argument('-t', '--tr_annot', required=True, help='TR annotation file')
	parser_tr.add_argument('-v', '--vcf_in', required=True, help='dummy input VCF file to use for header')
	parser_tr.add_argument('-o', '--vcf_out', required=True, help='output VCF file')
	parser_tr.add_argument('-s', '--sample_bam_file', required=True, help='a map between samples and bam files as a tab delimited text file. First column is the samples, and second column is the absolute path to the bam files')
	parser_tr.add_argument('-c', '--contig', default=None, help='contig name. If used the input VCF should be indexed')
	parser_tr.add_argument('-n', '--n_section', type=int, default=1, help='number of sections in the input VCF file')
	parser_tr.add_argument('-i', '--i_section', type=int, default=0, help='which section of the input VCF file to process')
	parser_tr.add_argument('-l', '--l_max', type=int, default=int(1e9), help='maximum TR span for genotyping')
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
