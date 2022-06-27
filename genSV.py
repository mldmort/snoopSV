import pysam
import glob
import pandas as pd
import numpy as np
import sys
from math import log10, factorial
import subprocess
from alignment.local_alignment import AlignmentScore

class sv_class:
	def __init__(self, svtype, chrom, start, stop, chr2, svlen):
		self.chrom = chrom
		self.start = start
		self.stop = stop
		self.info = {'SVTYPE':svtype, 'SVLEN':svlen, 'CHR2':chr2}

def infer_gt_sv(DR, DV, p_err):
	
	try:
		p_D_01 = (0.5**DR) * (0.5**DV)
	except:
		print('probability calculation problem p_D_01, DR:', DR, 'DV:', DV)
		p_D_01 = 0.0

	try:
		p_D_11 = (p_err**DR) * ((1.-p_err)**DV)
	except:
		print('probability calculation problem p_D_11, DR:', DR, 'DV:', DV)
		p_D_11 = 0.0

	try:
		p_D_00 = ((1.-p_err)**DR) * (p_err**DV)
	except:
		print('probability calculation problem p_D_00, DR:', DR, 'DV:', DV)
		p_D_00 = 0.0
	
	p_G = 1./3.
	p_D = (p_D_00 + p_D_01 + p_D_11)*p_G
	if p_D == 0:
		return './.', 0, 0, 0, 0, 0
		#print('p_D_00:', p_D_00)
		#print('p_D_01:', p_D_01)
		#print('p_D_11:', p_D_11)
		#print('DR:', DR)
		#print('DV:', DV)
	p_01_D = p_D_01 * p_G / p_D
	p_11_D = p_D_11 * p_G / p_D
	p_00_D = p_D_00 * p_G / p_D

	probs = [(p_01_D, '0/1', 1), (p_11_D, '1/1', 2), (p_00_D, '0/0', 0)]
	probs_sorted = sorted(probs, key=lambda x: x[0], reverse=True)
	GQ = round(-10.*log10(max(1. - probs_sorted[0][0], 1e-100)))
	SQ = round(-10.*log10(max(p_00_D, 1e-100)))

	return probs_sorted[0][1], GQ, p_11_D, p_01_D, p_00_D, SQ

def infer_gt_tr(count_list, p_err, svtype):
	if len(count_list) == 0:
		return './.', 0
	r = 2.*p_err/(1.-p_err)
	#print('r:', r)
	count_set = set(count_list)
	count_uniq_sorted_list = sorted(list(count_set))
	#print('count_uniq_sorted_list:', count_uniq_sorted_list)
	genotype_list = []
	for i in range(len(count_uniq_sorted_list)):
		for j in range(i,len(count_uniq_sorted_list)):
			genotype_list.append((count_uniq_sorted_list[i], count_uniq_sorted_list[j]))
	#print('genotype_list:', genotype_list)

	prob_C = {}
	for c_tuple in genotype_list:
		ci, cj = c_tuple
		prob_C[(ci,cj)] = 1.
		for c in count_list:
			if ci==cj and c==ci:
				prob_c = 1. - r
			elif c==ci:
				prob_c = 0.5*((1.-r)+p_err**(abs(c-cj)))
			elif c==cj:
				prob_c = 0.5*((1.-r)+p_err**(abs(c-ci)))
			else:
				prob_c = 0.5*(p_err**(abs(c-ci)) + p_err**(abs(c-cj)))
			prob_C[(ci,cj)] *= prob_c
	#print('prob_C:', prob_C)
	prob_C_sum = 0
	for c_tuple in genotype_list:
		ci, cj = c_tuple
		prob_C_sum += prob_C[(ci,cj)]
	if prob_C_sum==0:
		return './.', 0
	prob_gt = {}
	for c_tuple in genotype_list:
		ci, cj = c_tuple
		prob_gt[(ci,cj)] = prob_C[(ci,cj)]/prob_C_sum
	#print('prob_gt:', prob_gt)
	genotype, prob = sorted(prob_gt.items(), key=lambda x:x[1], reverse=True)[0]
	GT = str(genotype[0])+'/'+str(genotype[1])
	GQ = round(-10.*log10(max(1. - prob, 1e-100)))
	
	return GT, GQ

def infer_gt_tr_phased(count_dict, r_min):
	
	h1_list = []
	h2_list = []
	h0_list = []
	h1_list.extend(count_dict['H1'])
	h2_list.extend(count_dict['H2'])
	h0_list.extend(count_dict['H0'])
	#if len(h0_list) > 0:
	#	if len(h1_list) == 0:
	#		h1_list.extend(h0_list)
	#	if len(h2_list) == 0:
	#		h2_list.extend(h0_list)

	if len(h1_list) >= r_min:
		h1_gt = str(int(np.median(h1_list)))
	else:
		h1_gt = '.'

	if len(h2_list) >= r_min:
		h2_gt = str(int(np.median(h2_list)))
	else:
		h2_gt = '.'

	return h1_gt+'|'+h2_gt

def get_cigar_dict(read, ins_len_thr, del_len_thr):
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_al_start = read.query_alignment_start
	read_al_stop = read.query_alignment_end
	cigar_t = read.cigartuples
	ind_start = -1
	ind_end = -1
	for i, c in enumerate(cigar_t):
		if c[0] != 4:
			ind_start = i
			break
	for i in range(len(cigar_t)-1,-1,-1):
		c = cigar_t[i]
		if c[0] != 4:
			ind_end = i
			break
	assert ind_start != -1, 'wrong ind_start '+str(cigar_t)
	assert ind_end != -1, 'wrong ind_end '+str(cigar_t) 
	cigar_dict = {'I':{'ref_pos':[], 'len':[], 'seq_pos':[]}, 'D':{'ref_pos':[]}}
	cur_ref_pos = read_ref_start
	cur_seq_pos = 0
	for c in cigar_t[ind_start:ind_end+1]:
		if c[0] == 0:
			cur_ref_pos += c[1]
			cur_seq_pos += c[1]
		elif c[0] == 1:
			if c[1] >= ins_len_thr:
				cigar_dict['I']['ref_pos'].append(cur_ref_pos)
				cigar_dict['I']['seq_pos'].append((cur_seq_pos, cur_seq_pos+c[1]))
				cigar_dict['I']['len'].append(c[1])
			cur_seq_pos += c[1]
		elif c[0] == 2:
			if c[1] >= del_len_thr:
				cigar_dict['D']['ref_pos'].append((cur_ref_pos, cur_ref_pos+c[1]))
			cur_ref_pos += c[1]
	assert cur_ref_pos == read_ref_stop, 'something wrong with CIGAR length addition, cur_ref_pos: '+str(cur_ref_pos)+', read_ref_stop: '+str(read_ref_stop)
	assert cur_seq_pos == read_al_stop-read_al_start, 'something wrong with CIGAR length addition, cur_seq_pos: '+str(cur_seq_pos)+', read_al_stop-read_al_start: '+str(read_al_stop-read_al_start)

	return cigar_dict, ind_start, ind_end
	
def get_SA_cigar_dict(SA_cigar):
	SA_dict = {'S_left':0, 'S_right':0, 'M':0, 'D':0, 'I':0}
	c_cur = ''
	left_S_set = False
	for c in SA_cigar:
		if c == 'S':
			if not left_S_set:
				SA_dict['S_left'] = int(c_cur)
				left_S_set = True
			else:
				SA_dict['S_right'] = int(c_cur)
			c_cur = ''
		elif c == 'M':
			SA_dict['M'] = int(c_cur)
			c_cur = ''
			left_S_set = True
		elif c == 'D':
			SA_dict['D'] = int(c_cur)
			c_cur = ''
			left_S_set = True
		elif c == 'I':
			SA_dict['I'] = int(c_cur)
			c_cur = ''
			left_S_set = True
		else:
			c_cur += c

	return SA_dict

def get_seq_segment(read, ref_start, ref_stop):

	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_al_start = read.query_alignment_start ### this always starts from the left side of CIGAR, no mater + or - strand read
	read_al_stop = read.query_alignment_end
	read_al_len = read.query_alignment_length
	#print('read_ref_start:', read_ref_start)
	#print('read_ref_stop:', read_ref_stop)
	#print('read_al_start:', read_al_start)
	#print('read_al_stop:', read_al_stop)
	#print('read_al_len:', read_al_len)
	#print('ref_start:', ref_start)
	#print('ref_stop:', ref_stop)

	seq_start = -1
	seq_stop = -1
	blank_start = 0
	blank_stop = 0

	if (read_ref_start > ref_stop) or (read_ref_stop < ref_start):
		return '', 0, 0

	if read_ref_start > ref_start:
		seq_start = 0

	if read_ref_stop < ref_stop:
		seq_stop = read_al_len

	#print('seq_start:', seq_start)
	#print('seq_stop:', seq_stop)

	cigar_t = read.cigartuples
	#print('cigar_t:')
	#print(cigar_t)
	#print('len(cigar_t):', len(cigar_t))
	ind_start = -1
	ind_end = -1
	for i, c in enumerate(cigar_t):
		if c[0] != 4:
			ind_start = i
			break
	for i in range(len(cigar_t)-1,-1,-1):
		c = cigar_t[i]
		if c[0] != 4:
			ind_end = i
			break
	assert ind_start != -1, 'wrong ind_start '+str(cigar_t)
	assert ind_end != -1, 'wrong ind_end '+str(cigar_t) 
	#print('ind_start:', ind_start)
	#print('ind_end:', ind_end)
	#print('cigar_t[ind_start-1]', cigar_t[ind_start-1])
	#print('cigar_t[ind_start]', cigar_t[ind_start])
	#print('cigar_t[ind_end+1]', cigar_t[ind_end+1])

	### even for - strand reads the sequence in the bam file is reported as if it is on the + strand. So the CIGAR which starts from left to right is consistant with the sequence itself. so you don't need to reverse complement the sequences.
	cur_ref_pos = read_ref_start
	cur_seq_pos = 0
	for c in cigar_t[ind_start:ind_end+1]:
		if c[0] == 0: # Match
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'M'
		elif c[0] == 1: # Ins
			nxt_ref_pos = cur_ref_pos
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'I'
		elif c[0] == 2: # Del
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos
			cur_cigar = 'D'
		else:
			assert 0==1, 'wrong here, c[0]: '+str(c[0])
		#print('c:', c, 'cur_cigar:', cur_cigar, 'cur_seq_pos:', cur_seq_pos, 'nxt_seq_pos:', nxt_seq_pos, 'cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos)

		if (ref_start >= cur_ref_pos) and (ref_start <= nxt_ref_pos) and (seq_start == -1):
			if (cur_cigar == 'M'):
				#print('seq_start set at M, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_start:', ref_start)
				seq_start = cur_seq_pos + (ref_start - cur_ref_pos)
			elif (cur_cigar == 'D'):
				#print('seq_start set at D/I, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_start:', ref_start, 'cur_cigar:', cur_cigar)
				seq_start = cur_seq_pos
				blank_start = nxt_ref_pos - ref_start
			else:
				assert 0==1, 'wrong cur_cigar: '+cur_cigar+', cur_cigar should not be I. should be cought earlier'
		if (ref_stop >= cur_ref_pos) and (ref_stop <= nxt_ref_pos) and (seq_stop == -1):
			if (cur_cigar == 'M'):
				#print('seq_stop set at M, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_stop:', ref_stop)
				seq_stop = cur_seq_pos + (ref_stop - cur_ref_pos)
			elif (cur_cigar == 'D'):
				#print('seq_start set at D/I, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_stop:', ref_stop, 'cur_cigar:', cur_cigar)
				seq_stop = cur_seq_pos
				blank_stop = ref_stop - cur_ref_pos
			else:
				assert 0==1, 'wrong cur_cigar: '+cur_cigar+', cur_cigar should not be I. should be cought earlier'

		cur_ref_pos = nxt_ref_pos
		cur_seq_pos = nxt_seq_pos

	assert cur_ref_pos == read_ref_stop, 'something wrong with CIGAR length addition, cur_ref_pos: '+str(cur_ref_pos)+', read_ref_stop: '+str(read_ref_stop)
	assert cur_seq_pos == read_al_stop-read_al_start, 'something wrong with CIGAR length addition, cur_seq_pos: '+str(cur_seq_pos)+', read_al_stop-read_al_start: '+str(read_al_stop-read_al_start)

	assert ((seq_start != -1) and (seq_stop != -1)), 'problem with sequence index, seq_start: '+str(seq_start)+', seq_stop: '+str(seq_stop)

	#print('seq_start:', seq_start)
	#print('seq_stop:', seq_stop)

	### even for - strand reads the sequence in the bam file is reported as if it is on the + strand
	### query_alignment_sequence: just have the aligned sequence in the region not the soft clipped sequences
	### query_sequence: have all the read sequence, including the soft clipped sequence
	### read_al_len: the length of query_alignment_sequence
	sequence = read.query_alignment_sequence[max(seq_start,0):min(seq_stop,read_al_len)]

	return sequence, blank_start, blank_stop
	
def get_seq_segment_supp(read_a, read_b, ref_start, ref_stop):

	#print('ref_start:', ref_start)
	#print('ref_stop:', ref_stop)

	read_a_ref_start = read_a.reference_start
	read_a_ref_stop = read_a.reference_end
	read_a_al_start = read_a.query_alignment_start ### this always starts from the left side of CIGAR, no mater + or - strand read
	read_a_al_stop = read_a.query_alignment_end
	read_a_al_len = read_a.query_alignment_length
	#print('read_a_ref_start:', read_a_ref_start)
	#print('read_a_ref_stop:', read_a_ref_stop)
	#print('read_a_al_start:', read_a_al_start)
	#print('read_a_al_stop:', read_a_al_stop)
	#print('read_a_al_len:', read_a_al_len)
	#print('len(read_a.query_sequence)', len(read_a.query_sequence))

	read_b_ref_start = read_b.reference_start
	read_b_ref_stop = read_b.reference_end
	read_b_al_start = read_b.query_alignment_start ### this always starts from the left side of CIGAR, no mater + or - strand read
	read_b_al_stop = read_b.query_alignment_end
	read_b_al_len = read_b.query_alignment_length
	#print('read_b_ref_start:', read_b_ref_start)
	#print('read_b_ref_stop:', read_b_ref_stop)
	#print('read_b_al_start:', read_b_al_start)
	#print('read_b_al_stop:', read_b_al_stop)
	#print('read_b_al_len:', read_b_al_len)
	#print('len(read_b.query_sequence)', len(read_b.query_sequence))

	seq_start = -1
	seq_stop = -1
	blank_start = 0
	blank_stop = 0

	cigar_t_a = read_a.cigartuples
	ind_start_a = -1
	ind_end_a = -1
	for i, c in enumerate(cigar_t_a):
		if c[0] != 4:
			ind_start_a = i
			break
	for i in range(len(cigar_t_a)-1,-1,-1):
		c = cigar_t_a[i]
		if c[0] != 4:
			ind_end_a = i
			break
	assert ind_start_a != -1, 'wrong ind_start '+str(cigar_t_a)
	assert ind_end_a != -1, 'wrong ind_end '+str(cigar_t_a) 
	#print('ind_start_a:', ind_start_a)
	#print('nd_end_a:', ind_end_a)
	#print('cigar_t_a[ind_start-1]', cigar_t_a[ind_start_a-1])
	#print('cigar_t_a[ind_start]', cigar_t_a[ind_start_a])
	#print('cigar_t_a[ind_end+1]', cigar_t_a[ind_end_a+1])

	cigar_t_b = read_b.cigartuples
	ind_start_b = -1
	ind_end_b = -1
	for i, c in enumerate(cigar_t_b):
		if c[0] != 4:
			ind_start_b = i
			break
	for i in range(len(cigar_t_b)-1,-1,-1):
		c = cigar_t_b[i]
		if c[0] != 4:
			ind_end_b = i
			break
	assert ind_start_b != -1, 'wrong ind_start '+str(cigar_t_b)
	assert ind_end_b != -1, 'wrong ind_end '+str(cigar_t_b) 
	#print('ind_start_b:', ind_start_b)
	#print('nd_end_b:', ind_end_b)
	#print('cigar_t_b[ind_start-1]', cigar_t_b[ind_start_b-1])
	#print('cigar_t_b[ind_start]', cigar_t_b[ind_start_b])
	#print('cigar_t_b[ind_end+1]', cigar_t_b[ind_end_b+1])

	### even for - strand reads the sequence in the bam file is reported as if it is on the + strand. So the CIGAR which starts from left to right is consistant with the sequence itself. so you don't need to reverse complement the sequences.
	cur_ref_pos = read_a_ref_start
	cur_seq_pos = 0
	for c in cigar_t_a[ind_start_a:ind_end_a+1]:
		if c[0] == 0: # Match
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'M'
		elif c[0] == 1: # Ins
			nxt_ref_pos = cur_ref_pos
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'I'
		elif c[0] == 2: # Del
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos
			cur_cigar = 'D'
		else:
			assert 0==1, 'wrong here, c[0]: '+str(c[0])

		if (ref_start >= cur_ref_pos) and (ref_start <= nxt_ref_pos) and (seq_start == -1):
			if (cur_cigar == 'M'):
				#print('seq_start set at M, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_start:', ref_start)
				seq_start = cur_seq_pos + (ref_start - cur_ref_pos)
			elif (cur_cigar == 'D'):
				#print('seq_start set at D/I, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_start:', ref_start, 'cur_cigar:', cur_cigar)
				seq_start = cur_seq_pos
				blank_start = nxt_ref_pos - ref_start
			else:
				assert 0==1, 'wrong cur_cigar: '+cur_cigar+', cur_cigar should not be I. should be cought earlier'

		cur_ref_pos = nxt_ref_pos
		cur_seq_pos = nxt_seq_pos

	assert cur_ref_pos == read_a_ref_stop, 'something wrong with CIGAR length addition, cur_ref_pos: '+str(cur_ref_pos)+', read_a_ref_stop: '+str(read_a_ref_stop)
	assert cur_seq_pos == read_a_al_stop-read_a_al_start, 'something wrong with CIGAR length addition, cur_seq_pos: '+str(cur_seq_pos)+', read_a_al_stop-read_a_al_start: '+str(read_a_al_stop-read_a_al_start)

	assert (seq_start != -1), 'problem with sequence index, seq_start: '+str(seq_start)+', seq_stop: '+str(seq_stop)+', read: '+read_a.query_name

	cur_ref_pos = read_b_ref_start
	cur_seq_pos = 0
	for c in cigar_t_b[ind_start_b:ind_end_b+1]:
		if c[0] == 0: # Match
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'M'
		elif c[0] == 1: # Ins
			nxt_ref_pos = cur_ref_pos
			nxt_seq_pos = cur_seq_pos + c[1]
			cur_cigar = 'I'
		elif c[0] == 2: # Del
			nxt_ref_pos = cur_ref_pos + c[1]
			nxt_seq_pos = cur_seq_pos
			cur_cigar = 'D'
		else:
			assert 0==1, 'wrong here, c[0]: '+str(c[0])

		if (ref_stop >= cur_ref_pos) and (ref_stop <= nxt_ref_pos) and (seq_stop == -1):
			if (cur_cigar == 'M'):
				#print('seq_stop set at M, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_stop:', ref_stop)
				seq_stop = cur_seq_pos + (ref_stop - cur_ref_pos)
			elif (cur_cigar == 'D'):
				#print('seq_start set at D/I, cur_ref_pos:', cur_ref_pos, 'nxt_ref_pos:', nxt_ref_pos, 'ref_stop:', ref_stop, 'cur_cigar:', cur_cigar)
				seq_stop = cur_seq_pos
				blank_stop = ref_stop - cur_ref_pos
			else:
				assert 0==1, 'wrong cur_cigar: '+cur_cigar+', cur_cigar should not be I. should be cought earlier'

		cur_ref_pos = nxt_ref_pos
		cur_seq_pos = nxt_seq_pos

	assert cur_ref_pos == read_b_ref_stop, 'something wrong with CIGAR length addition, cur_ref_pos: '+str(cur_ref_pos)+', read_b_ref_stop: '+str(read_b_ref_stop)
	assert cur_seq_pos == read_b_al_stop-read_b_al_start, 'something wrong with CIGAR length addition, cur_seq_pos: '+str(cur_seq_pos)+', read_b_al_stop-read_b_al_start: '+str(read_b_al_stop-read_b_al_start)

	assert (seq_stop != -1), 'problem with sequence index, seq_start: '+str(seq_start)+', seq_stop: '+str(seq_stop)

	#print('seq_start:', seq_start)
	#print('seq_stop:', seq_stop)

	### even for - strand reads the sequence in the bam file is reported as if it is on the + strand
	### query_alignment_sequence: just have the aligned sequence in the region not the soft clipped sequences
	### query_sequence: have all the read sequence, including the soft clipped sequence
	### read_al_len: the length of query_alignment_sequence
	#sequence = read.query_alignment_sequence[max(seq_start,0):min(seq_stop,read_al_len)]
	sequence = read_a.query_sequence[(read_a_al_start+seq_start):(read_b_al_start+seq_stop)]

	return sequence, blank_start, blank_stop
	
def calc_recip_overlap(s1, e1, s2, e2):
	if (s1 > e2) or (s2 > e1):
		return 0
	so = max(s1, s2)
	eo = min(e1, e2)
	lo = eo - so
	assert lo >= 0, 'problem in calc_recip_overlap, s1,e1,s2,e2,lo: '+str(s1)+','+str(e1)+','+str(s2)+','+str(e2)+','+str(lo)
	f1 = float(lo) / float(e1 - s1)
	f2 = float(lo) / float(e2 - s2)
	f = min(f1, f2)
	return f

def sv_signiture(read, target_sv):

	CG_read_supp = False
	SA_read_supp = False
	locus_read = False

	CG_read_name = ''
	SA_read_name = ''
	locus_read_name = ''
	
	target_svtype = target_sv.info['SVTYPE']
	target_svlen = abs(target_sv.info['SVLEN'])
	target_start = target_sv.start
	target_stop = target_sv.stop
	
	### settings:
	mapping_quality_thr = 20
	d_max = 500
	region_buffer_left = target_start - d_max
	region_buffer_right = target_stop + d_max
	len_ratio_tol = 0.25
	ins_len_thr = 20
	del_len_thr = 20
	del_recip_overlap_thr = 0.3
	#ins_recip_overlap_thr = 0.01

	read_chrom = read.reference_name
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_ref_span = read_ref_stop - read_ref_start
	read_al_len = read.query_alignment_length
	read_al_start = read.query_alignment_start
	read_al_stop = read.query_alignment_end
	read_strand = '+'
	if read.is_reverse:
		read_strand = '-'

	#print('++++++++++++++++++++++++++++++')
	#print('target_start:', target_start)
	#print('target_stop:', target_stop)
	#print('read_ref_start:', read_ref_start)
	#print('read_ref_stop:', read_ref_stop)
	#print('read_strand:', read_strand)

	### 100 is an arbitrary length
	locus_region_len = 100
	if (target_start - locus_region_len > read_ref_start) and (target_start - locus_region_len < read_ref_stop):
		locus_read = True
		locus_read_name = read.query_name
	elif (target_start + locus_region_len > read_ref_start) and (target_start + locus_region_len < read_ref_stop):
		locus_read = True
		locus_read_name = read.query_name
	elif (target_stop - locus_region_len > read_ref_start) and (target_stop - locus_region_len < read_ref_stop):
		locus_read = True
		locus_read_name = read.query_name
	elif (target_stop + locus_region_len > read_ref_start) and (target_stop + locus_region_len < read_ref_stop):
		locus_read = True
		locus_read_name = read.query_name

	if not locus_read:
		return '', '', ''

	#print('read_ref_start:', read_ref_start)
	#print('read_ref_stop:', read_ref_stop)
	#print('read_ref_span:', read_ref_span)
	#print('read_al_len:', read_al_len)
	#print('read_strand:', read_strand)

	#print('read.query_alignment_length:')
	#print(read.query_alignment_length)
	#print('read.query_length:')
	#print(read.query_length)
	#print('len(read.query_sequence):')
	#print(len(read.query_sequence))
	#print('len(read.query_alignment_sequence):')
	#print(len(read.query_alignment_sequence))

	cigar_dict, ind_start, ind_end = get_cigar_dict(read, ins_len_thr, del_len_thr)

	#print('cigar_dict:', cigar_dict)
	#cigar_t = read.cigartuples
	#print('cigar_t[0]:', cigar_t[0])
	#print('cigar_t[-1]:', cigar_t[-1])

	##### IMPORTANT: if the read is on negative strand, read position starts from the left of CIGAR. CIGAR is always in the positive direction. pysam library doesn't respect this and always sets query_alignment_start to the left S value.
	cigar_t = read.cigartuples
	if read.is_reverse:
		if ind_end == len(cigar_t)-1:
			read_al_start = 0
		else:
			read_al_start = cigar_t[ind_end+1][1]
		read_al_stop = read_al_start + read_al_len

	### from CIGAR
	if target_svtype == 'INS':
		ref_pos_list = cigar_dict['I']['ref_pos']
		ref_len_list = cigar_dict['I']['len']
		sum_cigar_len = 0
		for ind, ref_pos in enumerate(ref_pos_list):
			if (ref_pos > region_buffer_left) and (ref_pos < region_buffer_right):
				#ref_len_half = int(ref_len_list[ind] / 2)
				#ref_pos_start = ref_pos - ref_len_half
				#ref_pos_stop = ref_pos + ref_len_half
				#recip_overlap = calc_recip_overlap(ref_pos_start, ref_pos_stop, target_start-ref_len_half, target_stop+ref_len_half)
				#if recip_overlap >= ins_recip_overlap_thr:
				#	CG_read_supp = True
				#	CG_read_name = read.query_name
				#	break
				ref_len = ref_len_list[ind]
				sum_cigar_len += ref_len
				if float(abs(ref_len-target_svlen))/float(target_svlen) < len_ratio_tol:
					CG_read_supp = True
					CG_read_name = read.query_name
					break
		#if (not CG_read_supp) and (float(abs(sum_cigar_len-target_svlen))/float(target_svlen) < len_ratio_tol):
		#	CG_read_supp = True
		#	CG_read_name = read.query_name
	elif target_svtype == 'DEL':
		ref_pos_list = cigar_dict['D']['ref_pos']
		sum_cigar_len = 0
		for ref_pos_t in ref_pos_list:
			ref_pos_start = ref_pos_t[0]
			ref_pos_stop = ref_pos_t[1]
			if (ref_pos_start > region_buffer_left) and (ref_pos_stop < region_buffer_right):
				recip_overlap = calc_recip_overlap(ref_pos_start, ref_pos_stop, target_start, target_stop)
				if recip_overlap >= del_recip_overlap_thr:
					CG_read_supp = True
					CG_read_name = read.query_name
					break
				#ref_len = ref_pos_stop - ref_pos_start
				#sum_cigar_len += ref_len
				#if float(abs(ref_len-target_svlen))/float(target_svlen) < len_ratio_tol:
				#	CG_read_supp = True
				#	CG_read_name = read.query_name
				#	break
		#if (not CG_read_supp) and (float(abs(sum_cigar_len-target_svlen))/float(target_svlen) < len_ratio_tol):
		#	CG_read_supp = True
		#	CG_read_name = read.query_name
	elif target_svtype == 'DUP':
		ref_pos_list = cigar_dict['I']['ref_pos']
		ref_len_list = cigar_dict['I']['len']
		for ind, ref_pos in enumerate(ref_pos_list):
			if (ref_pos > region_buffer_left) and (ref_pos < region_buffer_right):
				sv_len = ref_len_list[ind]
				if float(abs(sv_len-target_svlen))/float(target_svlen) < len_ratio_tol:
					CG_read_supp = True
					CG_read_name = read.query_name
					break

	#chr22,36701741,-,1433S4788M248D24S,60,725;	
	SA_next_right = {'SA_ref_start':-1, 'SA_ref_stop':-1, 'SA_read_start':1e15, 'SA_read_stop':1e15}
	SA_next_left = {'SA_ref_start':-1, 'SA_ref_stop':-1, 'SA_read_start':-1, 'SA_read_stop':-1}
	if read.has_tag("SA"):
		SA_tag = read.get_tag(tag="SA")
		SA_list = SA_tag.split(';')[:-1]
		for SA in SA_list:
			SA = SA.split(',')
			SA_chrom = SA[0]
			SA_ref_start = int(SA[1])
			SA_strand = SA[2]
			SA_cigar = SA[3]
			SA_mapq = float(SA[4])
			#print('SA_chrom:', SA_chrom)
			#print('SA_strand:', SA_strand)
			#print('SA_cigar:', SA_cigar)
			#print('SA_mapq:', SA_mapq)

			#### don't consider SA if:
			if (SA_mapq < mapping_quality_thr):
				continue
			if ((target_svtype == 'DEL') or (target_svtype == 'INS')) and (SA_strand != read_strand):
				continue
			if (target_svtype == 'TRA') and (SA_chrom == read_chrom):
				continue
			if (target_svtype != 'TRA') and (SA_chrom != read_chrom):
				continue
			if (target_svtype == 'INV') and (SA_strand == read_strand):
				continue
				
			SA_dict = {'S_left':0, 'S_right':0, 'M':0, 'D':0, 'I':0}
			c_cur = ''
			left_S_set = False
			for c in SA_cigar:
				if c == 'S':
					if not left_S_set:
						SA_dict['S_left'] = int(c_cur)
						left_S_set = True
					else:
						SA_dict['S_right'] = int(c_cur)
					c_cur = ''
				elif c == 'M':
					SA_dict['M'] = int(c_cur)
					c_cur = ''
					left_S_set = True
				elif c == 'D':
					SA_dict['D'] = int(c_cur)
					c_cur = ''
					left_S_set = True
				elif c == 'I':
					SA_dict['I'] = int(c_cur)
					c_cur = ''
					left_S_set = True
				else:
					c_cur += c

			SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']

			if (SA_strand == '+'):
				SA_read_start = SA_dict['S_left']
				SA_read_stop = SA_read_start + SA_dict['M'] + SA_dict['I']
			else:
				SA_read_start = SA_dict['S_right']
				SA_read_stop = SA_read_start + SA_dict['M'] + SA_dict['I']

			if (read_al_start < SA_read_start):
				delta_read = SA_read_start - read_al_stop
				if (read_strand=='+') and (SA_strand == '+'):
					delta_ref = SA_ref_start - read_ref_stop
					ref_overlap = -1*delta_ref
					bp1_overlap = SA_ref_start
					bp2_overlap = read_ref_stop
				elif (read_strand=='-') and (SA_strand == '-'):
					delta_ref = read_ref_start - SA_ref_stop
					ref_overlap = -1*delta_ref
					bp1_overlap = read_ref_start
					bp2_overlap = SA_ref_stop
				else:
					ref_overlap = min(read_ref_stop, SA_ref_stop) - max(read_ref_start, SA_ref_start)
					bp1_overlap = max(read_ref_start, SA_ref_start)
					bp2_overlap = min(read_ref_stop, SA_ref_stop)
			else:
				delta_read = read_al_start - SA_read_stop
				if (read_strand=='+') and (SA_strand == '+'):
					delta_ref = read_ref_start - SA_ref_stop
					ref_overlap = -1*delta_ref
					bp1_overlap = read_ref_start
					bp2_overlap = SA_ref_stop
				elif (read_strand=='-') and (SA_strand == '-'):
					delta_ref = SA_ref_start - read_ref_stop
					ref_overlap = -1*delta_ref
					bp1_overlap = SA_ref_start
					bp2_overlap = read_ref_stop
				else:
					ref_overlap = min(read_ref_stop, SA_ref_stop) - max(read_ref_start, SA_ref_start)
					bp1_overlap = max(read_ref_start, SA_ref_start)
					bp2_overlap = min(read_ref_stop, SA_ref_stop)
					
			#print('SA_dict:', SA_dict)
			#print('SA_ref_start', SA_ref_start)
			#print('SA_ref_stop', SA_ref_stop)
			#print('SA_read_start', SA_read_start)
			#print('SA_read_stop', SA_read_stop)
			#print('read_ref_start', read_ref_start)
			#print('read_ref_stop', read_ref_stop)
			#print('read_al_start', read_al_start)
			#print('read_al_stop', read_al_stop)
			#print('delta_ref:', delta_ref)
			#print('delta_read:', delta_read)
			#print('ref_overlap:', ref_overlap)

			if target_svtype == 'DEL':
				sv_len = delta_ref - delta_read
				if (ref_overlap < 30) and (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol):
					SA_read_supp = True
					SA_read_name = read.query_name
					break
			elif target_svtype == 'INS':
				sv_len = delta_read - delta_ref
				if (ref_overlap < 30) and (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol):
					SA_read_supp = True
					SA_read_name = read.query_name
					break
			elif target_svtype == 'DUP':
				if (read_strand == SA_strand):
					sv_len = ref_overlap
					#if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
					#	(   (float(abs(target_start - bp1_overlap))/float(target_svlen) < len_ratio_tol) or \
					#		(float(abs(target_stop - bp2_overlap))/float(target_svlen) < len_ratio_tol) ):
					if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
						(   (float(abs(target_start - bp1_overlap)) < 100) or \
							(float(abs(target_stop - bp2_overlap)) < 100) ):
						SA_read_supp = True
						SA_read_name = read.query_name
						break
			elif target_svtype == 'INVDUP':
				if (read_strand != SA_strand):
					sv_len = ref_overlap
					#if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
					#	(   (float(abs(target_start - bp1_overlap))/float(target_svlen) < len_ratio_tol) or \
					#		(float(abs(target_stop - bp2_overlap))/float(target_svlen) < len_ratio_tol) ):
					if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
						(   (float(abs(target_start - bp1_overlap)) < 100) or \
							(float(abs(target_stop - bp2_overlap)) < 100) ):
						SA_read_supp = True
						SA_read_name = read.query_name
						break
			elif target_svtype == 'INV':
				if (SA_read_start > read_al_start) and (SA_read_start < SA_next_right['SA_read_start']): # SA is right of read
					SA_next_right['SA_read_start'] = SA_read_start
					SA_next_right['SA_read_stop'] = SA_read_stop
					SA_next_right['SA_ref_start'] = SA_ref_start
					SA_next_right['SA_ref_stop'] = SA_ref_stop
				if (SA_read_start < read_al_start) and (SA_read_start > SA_next_left['SA_read_start']): # SA is left of read
					SA_next_left['SA_read_start'] = SA_read_start
					SA_next_left['SA_read_stop'] = SA_read_stop
					SA_next_left['SA_ref_start'] = SA_ref_start
					SA_next_left['SA_ref_stop'] = SA_ref_stop
			elif target_svtype == 'TRA':
				if (read_chrom == target_sv.chrom) and \
					(SA_chrom == target_sv.info['CHR2']) and \
					( (abs(read_ref_stop - target_start) < 100) or (abs(read_ref_start - target_start) < 100) ) and \
					( (abs(SA_ref_stop - target_stop) < 100) or (abs(SA_ref_start - target_stop) < 100) ):
					SA_read_supp = True
					SA_read_name = read.query_name
					break
						
	if target_svtype == 'INV':
		#print('SA_next_right:', SA_next_right)
		#print('SA_next_left:', SA_next_left)
		### read has a right SA
		if (SA_next_right['SA_ref_start'] != -1):
			if read_strand == '+':
				breakpoint1 = read_ref_stop
				breakpoint2 = SA_next_right['SA_ref_stop']
			else:
				breakpoint1 = read_ref_start
				breakpoint2 = SA_next_right['SA_ref_start']
			min_bp_right = min(breakpoint1, breakpoint2)
			max_bp_right = max(breakpoint1, breakpoint2)
		else:
			min_bp_right = -1
			max_bp_right = 1e15
		### read has a left SA
		if (SA_next_left['SA_ref_start'] != -1):
			if read_strand == '+':
				breakpoint1 = read_ref_start
				breakpoint2 = SA_next_left['SA_ref_start']
			else:
				breakpoint1 = read_ref_stop
				breakpoint2 = SA_next_left['SA_ref_stop']
			min_bp_left = min(breakpoint1, breakpoint2)
			max_bp_left = max(breakpoint1, breakpoint2)
		else:
			min_bp_left = -1
			max_bp_left = 1e15
		
		bp_1 = max(min_bp_right, min_bp_left)
		bp_2 = min(max_bp_right, max_bp_left)
		#print('bp_1:', bp_1, 'bp_2:', bp_2)
		if bp_1 != -1:
			sv_len = abs(bp_2 - bp_1)
			if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol):
				SA_read_supp = True
				SA_read_name = read.query_name
		
	#print ('locus_read_name:', locus_read_name, 'CG_read_name:', CG_read_name, 'SA_read_name:', SA_read_name)
	return locus_read_name, CG_read_name, SA_read_name

def tr_signature(read, target_sv, tr_start_list, tr_end_list, period_len_list, CN_list, period_seq_list, k_s_dict, visited_read_set_list):
	
	#print('+++++++++++ str signature ++++++++++')
	#print('period_len_list:', period_len_list)
	#print('CN_list:', CN_list)
	#print('period_seq_list:', period_seq_list)
	#
	#print('read.query_alignment_start:', read.query_alignment_start)
	#print('read.query_alignment_end:', read.query_alignment_end)
	#print('read.query_alignment_length:', read.query_alignment_length)
	#print('read.query_length:', read.query_length)
	#print('len(read.query_sequence):', len(read.query_sequence))
	#print('len(read.query_alignment_sequence):', len(read.query_alignment_sequence))
	#print('read.reference_start:', read.reference_start)
	#print('read.reference_end:', read.reference_end)

	### list of read name and repeat number for each TR in TR list
	locus_read_name_list = []
	num_repeat_align_list = []
	num_repeat_length_list = []
	num_bp_list = []
	
	read_chrom = read.reference_name
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end

	for i_tr in range(len(tr_start_list)):

		tr_start = tr_start_list[i_tr]
		tr_end = tr_end_list[i_tr]
		period_len = period_len_list[i_tr]
		period_seq = period_seq_list[i_tr]
		flanking_bp = period_len
		visited_read_set = visited_read_set_list[i_tr]
		#print('tr_start:', tr_start)
		#print('tr_end:', tr_end)
		#print('period_len:', period_len)
		#print('period_seq:', period_seq)

		if ((read_ref_start <= tr_start) and (read_ref_stop >= tr_end) and (read.query_name not in visited_read_set)):

			locus_read_name_list.append(read.query_name)

			read_seq_tr, _, _ = get_seq_segment(read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
			#score_ind_list, raw_score_ind_list = AlignmentScore(period_seq, read_seq_tr, k_s_dict)
			#num_repeat_align = len(score_ind_list)
			num_repeat_align = -1
			num_repeat_length = int(round((len(read_seq_tr)-2.*flanking_bp)/period_len))
			num_repeat_align_list.append(num_repeat_align)
			num_repeat_length_list.append(num_repeat_length)
			num_bp_list.append(len(read_seq_tr)-2*flanking_bp)
			#print(read_seq_tr, num_repeat_align, num_repeat_length, read.query_name)

		else:
			locus_read_name_list.append('')
			num_repeat_align_list.append(-1)
			num_repeat_length_list.append(-1)
			num_bp_list.append(-1)

	return locus_read_name_list, num_repeat_align_list, num_repeat_length_list, num_bp_list

def tr_signature_2(read, tr_start, tr_end, period_len, CN, period_seq, k_s_dict, visited_read_set, bam_file, mapping_quality_thr):
	
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_strand = '+'
	if read.is_reverse:
		read_strand = '-'

	#flanking_bp = period_len
	flanking_bp = 20

	#if ((read_ref_start <= tr_start) and (read_ref_stop >= tr_end) and (read.query_name not in visited_read_set)):
	if (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):

		#print('read.query_name', read.query_name)
		#print('tr_start:', tr_start)
		#print('tr_end:', tr_end)
		#print('period_len:', period_len)
		#print('period_seq:', period_seq)

		locus_read_name = read.query_name

		read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		#read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start, tr_end) ### tr_start and stop are 0-based
		#print('len(read_seq_tr):')
		#print(len(read_seq_tr))
		#score_ind_list, raw_score_ind_list = AlignmentScore(period_seq, read_seq_tr, k_s_dict)
		#num_repeat_align = len(score_ind_list)
		num_repeat_align = -1
		#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
		num_repeat_length = -1
		num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
		#num_bp = len(read_seq_tr) - blank_start - blank_stop
		#print(read_seq_tr, num_repeat_align, num_repeat_length, read.query_name)

	elif (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_start - flanking_bp)) and (read.query_name not in visited_read_set):
		SA_reads = []
		has_sa_read = False
		if read.has_tag("SA"):
			SA_tag = read.get_tag(tag="SA")
			SA_list = SA_tag.split(';')[:-1]
			for SA in SA_list:
				SA = SA.split(',')
				SA_chrom = SA[0]
				SA_ref_start = int(SA[1])
				SA_strand = SA[2]
				SA_cigar = SA[3]
				SA_mapq = float(SA[4])
				SA_dict = get_SA_cigar_dict(SA_cigar)
				SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']
				if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_end + flanking_bp)) and (SA_ref_stop > (tr_end + flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
					has_sa_read = True
					break
		if has_sa_read:
			fh_bam = pysam.AlignmentFile(bam_file, 'rb')
			for sa_read in fh_bam.fetch(read.reference_name, max(0,tr_start-flanking_bp), tr_end+flanking_bp):
				sa_read_ref_start = sa_read.reference_start
				sa_read_ref_stop = sa_read.reference_end
				sa_read_strand = '+'
				if sa_read.is_reverse:
					sa_read_strand = '-'
				sa_read_mapq = sa_read.mapping_quality
				if (sa_read.query_name == read.query_name) and (sa_read_ref_start < (tr_end + flanking_bp)) and (sa_read_ref_stop > (tr_end + flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr):
					SA_reads.append(sa_read)

		locus_read_name = read.query_name

		if len(SA_reads)>0:
			sa_read = SA_reads[0]

			read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
			#read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start, tr_end) ### tr_start and stop are 0-based
			num_repeat_align = -1
			#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
			num_repeat_length = -1
			num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
			#num_bp = len(read_seq_tr) - blank_start - blank_stop
		else:
			num_repeat_align = -1
			num_repeat_length = -1
			num_bp = -1
	elif (read_ref_start < (tr_end + flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):
		SA_reads = []
		has_sa_read = False
		if read.has_tag("SA"):
			SA_tag = read.get_tag(tag="SA")
			SA_list = SA_tag.split(';')[:-1]
			for SA in SA_list:
				SA = SA.split(',')
				SA_chrom = SA[0]
				SA_ref_start = int(SA[1])
				SA_strand = SA[2]
				SA_cigar = SA[3]
				SA_mapq = float(SA[4])
				SA_dict = get_SA_cigar_dict(SA_cigar)
				SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']
				if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_start - flanking_bp)) and (SA_ref_stop > (tr_start - flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
					has_sa_read = True
					break
		if has_sa_read:
			fh_bam = pysam.AlignmentFile(bam_file, 'rb')
			for sa_read in fh_bam.fetch(read.reference_name, max(0,tr_start-flanking_bp), tr_end+flanking_bp):
				sa_read_ref_start = sa_read.reference_start
				sa_read_ref_stop = sa_read.reference_end
				sa_read_strand = '+'
				if sa_read.is_reverse:
					sa_read_strand = '-'
				sa_read_mapq = sa_read.mapping_quality
				if (sa_read.query_name == read.query_name) and (sa_read_ref_start < (tr_start - flanking_bp)) and (sa_read_ref_stop > (tr_start - flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr):
					SA_reads.append(sa_read)

		locus_read_name = read.query_name

		if len(SA_reads)>0:
			sa_read = SA_reads[0]

			read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(sa_read, read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
			#read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(sa_read, read, tr_start, tr_end) ### tr_start and stop are 0-based
			num_repeat_align = -1
			#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
			num_repeat_length = -1
			num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
			#num_bp = len(read_seq_tr) - blank_start - blank_stop
		else:
			num_repeat_align = -1
			num_repeat_length = -1
			num_bp = -1
	else:
		locus_read_name = ''
		num_repeat_align = -1
		num_repeat_length = -1
		num_bp = -1

	return locus_read_name, num_repeat_align, num_repeat_length, num_bp

def tr_signature_3(read, tr_start, tr_end, period_len, CN, period_seq, k_s_dict, visited_read_set, bam_file, mapping_quality_thr):
	
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_strand = '+'
	if read.is_reverse:
		read_strand = '-'

	#flanking_bp = period_len
	flanking_bp = 20

	#print('read.query_name', read.query_name)
	#print('tr_start:', tr_start)
	#print('tr_end:', tr_end)
	#print('read_ref_start:', read_ref_start)
	#print('read_ref_stop:', read_ref_stop)
	#print('read.query_alignment_start:', read.query_alignment_start)
	#print('read.query_alignment_end:', read.query_alignment_end)

	#### find SA alignments which span tr_end+flanking_bp. We check if read spans tr_start-flanking_bp later
	#### we assume pysam fetch pulls reads from smallest ref position
	has_sa_read = False
	if read.has_tag("SA"):
		SA_tag = read.get_tag(tag="SA")
		SA_list = SA_tag.split(';')[:-1]
		for SA in SA_list:
			SA = SA.split(',')
			SA_chrom = SA[0]
			SA_ref_start = int(SA[1]) - 1 # bam is 1-based, we turn it into 0-based to be consistant with pysam
			SA_strand = SA[2]
			SA_cigar = SA[3]
			SA_mapq = float(SA[4])
			SA_dict = get_SA_cigar_dict(SA_cigar)
			SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']
			#print('SA:', SA)
			#print('SA_chrom:', SA_chrom)
			#print('SA_ref_start:', SA_ref_start)
			#print('SA_ref_stop:', SA_ref_stop)
			#print('SA_strand:', SA_strand)
			#print('SA_cigar:', SA_cigar)
			#print('SA_mapq:', SA_mapq)
			if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_end + flanking_bp)) and (SA_ref_stop > (tr_end + flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
				has_sa_read = True
				break

	SA_reads = []
	if has_sa_read:
		key_read_self = '_'.join([str(read_ref_start), str(read_ref_stop), read.get_tag(tag="SA")])
		fh_bam = pysam.AlignmentFile(bam_file, 'rb')
		for sa_read in fh_bam.fetch(read.reference_name, max(0,tr_start-flanking_bp), tr_end+flanking_bp):
			if (sa_read.query_name == read.query_name):
				sa_read_ref_start = sa_read.reference_start
				sa_read_ref_stop = sa_read.reference_end
				sa_read_strand = '+'
				if sa_read.is_reverse:
					sa_read_strand = '-'
				sa_read_mapq = sa_read.mapping_quality
				sa_read_SA = ''
				if sa_read.has_tag("SA"):
					sa_read_SA = sa_read.get_tag(tag="SA")
				key_sa_read = '_'.join([str(sa_read_ref_start), str(sa_read_ref_stop), sa_read_SA])
				if (sa_read_ref_start < (tr_end + flanking_bp)) and (sa_read_ref_stop > (tr_end + flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr) and (key_sa_read != key_read_self):
					SA_reads.append(sa_read)
					#print('sa_read_ref_start:', sa_read_ref_start)
					#print('sa_read_ref_stop:', sa_read_ref_stop)
					#print('sa_read_strand:', sa_read_strand)
					#print('sa_read_mapq:', sa_read_mapq)
					#print('sa_read.query_alignment_start:', sa_read.query_alignment_start)
					#print('sa_read.query_alignment_end:', sa_read.query_alignment_end)

	### query_alignment_start: this always starts from the left side of CIGAR, no mater + or - strand read
	SA_reads.sort(key=lambda x: x.query_alignment_start) # sorts from small to large, the last one is the largest
	#print('SA_reads:', SA_reads)

	if (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_start - flanking_bp)) and (len(SA_reads)>0) and (read.query_name not in visited_read_set):

		locus_read_name = read.query_name
		sa_read = SA_reads[-1]

		#read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(sa_read, read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		num_repeat_align = -1
		#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
		num_repeat_length = -1
		num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
	elif (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):

		#print('read.query_name', read.query_name)
		#print('tr_start:', tr_start)
		#print('tr_end:', tr_end)
		#print('period_len:', period_len)
		#print('period_seq:', period_seq)

		locus_read_name = read.query_name

		read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		#print('len(read_seq_tr):')
		#print(len(read_seq_tr))
		#score_ind_list, raw_score_ind_list = AlignmentScore(period_seq, read_seq_tr, k_s_dict)
		#num_repeat_align = len(score_ind_list)
		num_repeat_align = -1
		#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
		num_repeat_length = -1
		num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
		#print(read_seq_tr, num_repeat_align, num_repeat_length, read.query_name)
	else:
		locus_read_name = ''
		num_repeat_align = -1
		num_repeat_length = -1
		num_bp = -1





	#**#if ((read_ref_start <= tr_start) and (read_ref_stop >= tr_end) and (read.query_name not in visited_read_set)):
	#**if (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):

	#**	#print('read.query_name', read.query_name)
	#**	#print('tr_start:', tr_start)
	#**	#print('tr_end:', tr_end)
	#**	#print('period_len:', period_len)
	#**	#print('period_seq:', period_seq)

	#**	locus_read_name = read.query_name

	#**	read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
	#**	#read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start, tr_end) ### tr_start and stop are 0-based
	#**	#print('len(read_seq_tr):')
	#**	#print(len(read_seq_tr))
	#**	#score_ind_list, raw_score_ind_list = AlignmentScore(period_seq, read_seq_tr, k_s_dict)
	#**	#num_repeat_align = len(score_ind_list)
	#**	num_repeat_align = -1
	#**	#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
	#**	num_repeat_length = -1
	#**	num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
	#**	#num_bp = len(read_seq_tr) - blank_start - blank_stop
	#**	#print(read_seq_tr, num_repeat_align, num_repeat_length, read.query_name)

	#**elif (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_start - flanking_bp)) and (read.query_name not in visited_read_set):
	#**	SA_reads = []
	#**	has_sa_read = False
	#**	if read.has_tag("SA"):
	#**		SA_tag = read.get_tag(tag="SA")
	#**		SA_list = SA_tag.split(';')[:-1]
	#**		for SA in SA_list:
	#**			SA = SA.split(',')
	#**			SA_chrom = SA[0]
	#**			SA_ref_start = int(SA[1])
	#**			SA_strand = SA[2]
	#**			SA_cigar = SA[3]
	#**			SA_mapq = float(SA[4])
	#**			SA_dict = get_SA_cigar_dict(SA_cigar)
	#**			SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']
	#**			if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_end + flanking_bp)) and (SA_ref_stop > (tr_end + flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
	#**				has_sa_read = True
	#**				break
	#**	if has_sa_read:
	#**		fh_bam = pysam.AlignmentFile(bam_file, 'rb')
	#**		for sa_read in fh_bam.fetch(read.reference_name, max(0,tr_start-flanking_bp), tr_end+flanking_bp):
	#**			sa_read_ref_start = sa_read.reference_start
	#**			sa_read_ref_stop = sa_read.reference_end
	#**			sa_read_strand = '+'
	#**			if sa_read.is_reverse:
	#**				sa_read_strand = '-'
	#**			sa_read_mapq = sa_read.mapping_quality
	#**			if (sa_read.query_name == read.query_name) and (sa_read_ref_start < (tr_end + flanking_bp)) and (sa_read_ref_stop > (tr_end + flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr):
	#**				SA_reads.append(sa_read)

	#**	locus_read_name = read.query_name

	#**	if len(SA_reads)>0:
	#**		sa_read = SA_reads[0]

	#**		read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
	#**		#read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start, tr_end) ### tr_start and stop are 0-based
	#**		num_repeat_align = -1
	#**		#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
	#**		num_repeat_length = -1
	#**		num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
	#**		#num_bp = len(read_seq_tr) - blank_start - blank_stop
	#**	else:
	#**		num_repeat_align = -1
	#**		num_repeat_length = -1
	#**		num_bp = -1
	#**elif (read_ref_start < (tr_end + flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):
	#**	SA_reads = []
	#**	has_sa_read = False
	#**	if read.has_tag("SA"):
	#**		SA_tag = read.get_tag(tag="SA")
	#**		SA_list = SA_tag.split(';')[:-1]
	#**		for SA in SA_list:
	#**			SA = SA.split(',')
	#**			SA_chrom = SA[0]
	#**			SA_ref_start = int(SA[1])
	#**			SA_strand = SA[2]
	#**			SA_cigar = SA[3]
	#**			SA_mapq = float(SA[4])
	#**			SA_dict = get_SA_cigar_dict(SA_cigar)
	#**			SA_ref_stop = SA_ref_start + SA_dict['M'] + SA_dict['D']
	#**			if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_start - flanking_bp)) and (SA_ref_stop > (tr_start - flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
	#**				has_sa_read = True
	#**				break
	#**	if has_sa_read:
	#**		fh_bam = pysam.AlignmentFile(bam_file, 'rb')
	#**		for sa_read in fh_bam.fetch(read.reference_name, max(0,tr_start-flanking_bp), tr_end+flanking_bp):
	#**			sa_read_ref_start = sa_read.reference_start
	#**			sa_read_ref_stop = sa_read.reference_end
	#**			sa_read_strand = '+'
	#**			if sa_read.is_reverse:
	#**				sa_read_strand = '-'
	#**			sa_read_mapq = sa_read.mapping_quality
	#**			if (sa_read.query_name == read.query_name) and (sa_read_ref_start < (tr_start - flanking_bp)) and (sa_read_ref_stop > (tr_start - flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr):
	#**				SA_reads.append(sa_read)

	#**	locus_read_name = read.query_name

	#**	if len(SA_reads)>0:
	#**		sa_read = SA_reads[0]

	#**		read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(sa_read, read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
	#**		#read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(sa_read, read, tr_start, tr_end) ### tr_start and stop are 0-based
	#**		num_repeat_align = -1
	#**		#num_repeat_length = int(round((len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop))/float(period_len)))
	#**		num_repeat_length = -1
	#**		num_bp = len(read_seq_tr)-max(0,flanking_bp-blank_start)-max(0,flanking_bp-blank_stop)
	#**		#num_bp = len(read_seq_tr) - blank_start - blank_stop
	#**	else:
	#**		num_repeat_align = -1
	#**		num_repeat_length = -1
	#**		num_bp = -1
	#**else:
	#**	locus_read_name = ''
	#**	num_repeat_align = -1
	#**	num_repeat_length = -1
	#**	num_bp = -1

	return locus_read_name, num_repeat_align, num_repeat_length, num_bp
