from math import log10
from snoopsv.utils import get_cigar_dict, calc_recip_overlap

class sv_class:
	def __init__(self, rec):
		self.id = rec.id
		self.chrom = rec.chrom
		self.start = rec.start
		self.stop = rec.stop
		svtype = rec.info['SVTYPE']
		svlen = rec.info['SVLEN']
		if svtype == 'TRA':
			chr2 = rec.info['CHR2']
		else:
			chr2 = self.chrom
		self.info = {'SVTYPE': svtype, 'SVLEN': svlen, 'CHR2':chr2}

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

def sv_signature(read, target_sv):

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
				if float(abs(ref_len-target_svlen))/float(target_svlen) < len_ratio_tol:
					CG_read_supp = True
					CG_read_name = read.query_name
					break
	elif target_svtype == 'DEL':
		ref_pos_list = cigar_dict['D']['ref_pos']
		for ref_pos_t in ref_pos_list:
			ref_pos_start = ref_pos_t[0]
			ref_pos_stop = ref_pos_t[1]
			if (ref_pos_start > region_buffer_left) and (ref_pos_stop < region_buffer_right):
				recip_overlap = calc_recip_overlap(ref_pos_start, ref_pos_stop, target_start, target_stop)
				if recip_overlap >= del_recip_overlap_thr:
					CG_read_supp = True
					CG_read_name = read.query_name
					break
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
				# check if the read has a DUP signal
				sv_len = ref_overlap
				if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
					(   (float(abs(target_start - bp1_overlap)) < 100) or \
						(float(abs(target_stop - bp2_overlap)) < 100) ):
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
