from math import log10
from snoopsv.utils import get_cigar_dict, calc_recip_overlap

class sv_class:

	def __init__(self, rec, len_ratio_tol, ins_len_thr, del_len_thr, del_recip_overlap_thr):

		self.len_ratio_tol = len_ratio_tol
		self.ins_len_thr = ins_len_thr
		self.del_len_thr = del_len_thr
		self.del_recip_overlap_thr = del_recip_overlap_thr

		self.id = rec.id
		self.chrom = rec.chrom
		self.start = rec.start
		self.stop = rec.stop

		try:
			self.svtype = rec.info['SVTYPE']
		except ValueError:
			raise ValueError(f'SVTYPE is not defined for id: {rec.id}. Check the VCF file!')

		if self.svtype not in ('BND', 'TRA'):
			try:
				self.svlen = abs(rec.info['SVLEN'])
			except ValueError:
				raise ValueError(f'SVLEN is not defined for id: {rec.id}. Check the VCF file!')
		else:
			self.svlen = None

		if self.svtype in ('BND', 'TRA'):
			try:
				self.chr2 = rec.info['CHR2']
			except ValueError:
				raise ValueError(f'CHR2 is not defined for id: {rec.id}. Check the VCF file!')
		else:
			self.chr2 = self.chrom

		if self.svtype == 'BND':
			if '[' in rec.alts[0]:
				chr2_pos2 = rec.alts[0].split('[')[1] # always second member, both N[... and [...
				chr2_pos2 = chr2_pos2.split(':')
				assert chr2_pos2[0] == self.chr2, f'problem with parsing alt in BND call, sv id: {sv_id}'
				self.pos2 = chr2_pos2[1]
			elif ']' in rec.alts[0]:
				chr2_pos2 = rec.alts[0].split(']')[1] # always second member, both N]... and ]...
				chr2_pos2 = chr2_pos2.split(':')
				assert chr2_pos2[0] == self.chr2, f'problem with parsing alt in BND call, sv id: {sv_id}'
				self.pos2 = chr2_pos2[1]
			else:
				raise NameError(f'[ or ] not found in alt column of BND svtype, for sv id: {sv_id}')
		else:
			self.pos2 = None

	def sv_len_pass(self, proposed_len, target_svlen):
		if abs(proposed_len - target_svlen) / target_svlen < self.len_ratio_tol:
			return True
		return False

	def set_read_variables(self, read, buffer_length):
		self.CG_read_supp = False
		self.SA_read_supp = False
		self.locus_read = False

		self.CG_read_name = ''
		self.SA_read_name = ''
		self.locus_read_name = ''

		self.region_buffer_left = self.start - buffer_length
		self.region_buffer_right = self.stop + buffer_length

		self.read_chrom = read.reference_name
		self.read_ref_start = read.reference_start
		self.read_ref_stop = read.reference_end
		self.read_ref_span = self.read_ref_stop - self.read_ref_start
		self.read_al_len = read.query_alignment_length
		self.read_al_start = read.query_alignment_start
		self.read_al_stop = read.query_alignment_end
		self.read_strand = '+'
		if read.is_reverse:
			self.read_strand = '-'

	def set_cigar_variables(self, read):
		self.cigar_dict, ind_start, ind_end = get_cigar_dict(read, self.ins_len_thr, self.del_len_thr)

		##### IMPORTANT: if the read is on negative strand, read position starts from the left of CIGAR. CIGAR is always in the positive direction. pysam library doesn't respect this and always sets query_alignment_start to the left S value.
		# Now we want to go to the read coordinates, not the ref coordinates for read_al_start and read_al_end.
		# To be able to compare read segments with +/- strands on the same read coordinate.
		cigar_t = read.cigartuples
		if read.is_reverse:
			if ind_end == len(cigar_t)-1:
				self.read_al_start = 0
			else:
				self.read_al_start = cigar_t[ind_end+1][1]
			self.read_al_stop = self.read_al_start + self.read_al_len

	# this is a generator function going over the SA tag variables
	def get_sa_variables(self, read, mapping_quality_thr):
		#print(f'read name: {read.query_name}, sv type: {self.svtype}, sv start: {self.start}')
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
				if ((self.svtype == 'DEL') or (self.svtype == 'INS')) and (SA_strand != self.read_strand):
					continue
				if (self.svtype == 'TRA') and (SA_chrom == self.read_chrom):
					continue
				if (self.svtype != 'TRA') and (SA_chrom != self.read_chrom):
					continue
				if (self.svtype == 'INV') and (SA_strand == self.read_strand):
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

				# only INS and DEL svtypes use this variable, but be careful
				# if other svtypes want to use it, it is not set for + - or - + strands
				delta_ref = None
				if (self.read_al_start < SA_read_start):
					delta_read = SA_read_start - self.read_al_stop
					if (self.read_strand == '+') and (SA_strand == '+'):
						delta_ref = SA_ref_start - self.read_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = SA_ref_start
						bp2_overlap = self.read_ref_stop
					elif (self.read_strand == '-') and (SA_strand == '-'):
						delta_ref = self.read_ref_start - SA_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = self.read_ref_start
						bp2_overlap = SA_ref_stop
					else:
						ref_overlap = min(self.read_ref_stop, SA_ref_stop) - max(self.read_ref_start, SA_ref_start)
						bp1_overlap = max(self.read_ref_start, SA_ref_start)
						bp2_overlap = min(self.read_ref_stop, SA_ref_stop)
				else:
					delta_read = self.read_al_start - SA_read_stop
					if (self.read_strand == '+') and (SA_strand == '+'):
						delta_ref = self.read_ref_start - SA_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = self.read_ref_start
						bp2_overlap = SA_ref_stop
					elif (self.read_strand == '-') and (SA_strand == '-'):
						delta_ref = SA_ref_start - self.read_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = SA_ref_start
						bp2_overlap = self.read_ref_stop
					else:
						ref_overlap = min(self.read_ref_stop, SA_ref_stop) - max(self.read_ref_start, SA_ref_start)
						bp1_overlap = max(self.read_ref_start, SA_ref_start)
						bp2_overlap = min(self.read_ref_stop, SA_ref_stop)

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

				if (SA_read_start > self.read_al_start) and (SA_read_start < SA_next_right['SA_read_start']): # SA is right of read
					SA_next_right['SA_read_start'] = SA_read_start
					SA_next_right['SA_read_stop'] = SA_read_stop
					SA_next_right['SA_ref_start'] = SA_ref_start
					SA_next_right['SA_ref_stop'] = SA_ref_stop
				if (SA_read_start < self.read_al_start) and (SA_read_start > SA_next_left['SA_read_start']): # SA is left of read
					SA_next_left['SA_read_start'] = SA_read_start
					SA_next_left['SA_read_stop'] = SA_read_stop
					SA_next_left['SA_ref_start'] = SA_ref_start
					SA_next_left['SA_ref_stop'] = SA_ref_stop

				yield SA_strand, delta_read, delta_ref, ref_overlap, bp1_overlap, bp2_overlap, SA_next_left, SA_next_right
		else:
			yield tuple([None] * 8)

	def ins_signature(self, read, mapping_quality_thr):
		# from CIGAR
		pos_list = self.cigar_dict['I']['ref_pos']
		len_list = self.cigar_dict['I']['len']
		for ind, ref_pos in enumerate(pos_list):
			if (ref_pos > self.region_buffer_left) and (ref_pos < self.region_buffer_right):
				proposed_len = len_list[ind]
				if self.sv_len_pass(proposed_len, self.svlen):
					self.CG_read_supp = True
					self.CG_read_name = read.query_name
					break

		# from supplementary alignments
		for SA_strand, delta_read, delta_ref, ref_overlap, bp1_overlap, bp2_overlap, SA_next_left, SA_next_right in self.get_sa_variables(read, mapping_quality_thr):
			if SA_strand == None:
				break
			proposed_len = delta_read - delta_ref
			if (ref_overlap < 30) and self.sv_len_pass(proposed_len, self.svlen):
				self.SA_read_supp = True
				self.SA_read_name = read.query_name
				break
			# check if the read has a DUP signal
			proposed_len = ref_overlap
			if self.sv_len_pass(proposed_len, self.svlen) and \
				(   (float(abs(self.start - bp1_overlap)) < 100) or \
					(float(abs(self.stop - bp2_overlap)) < 100) ):
				self.SA_read_supp = True
				self.SA_read_name = read.query_name
				break

	def del_signature(self, read, mapping_quality_thr):
		# from CIGAR
		ref_pos_list = self.cigar_dict['D']['ref_pos']
		for ref_pos_t in ref_pos_list:
			ref_pos_start = ref_pos_t[0]
			ref_pos_stop = ref_pos_t[1]
			if (ref_pos_start > self.region_buffer_left) and (ref_pos_stop < self.region_buffer_right):
				recip_overlap = calc_recip_overlap(ref_pos_start, ref_pos_stop, self.start, self.stop)
				if recip_overlap >= self.del_recip_overlap_thr:
					self.CG_read_supp = True
					self.CG_read_name = read.query_name
					break

		# from supplementary alignments
		for SA_strand, delta_read, delta_ref, ref_overlap, bp1_overlap, bp2_overlap, SA_next_left, SA_next_right in self.get_sa_variables(read, mapping_quality_thr):
			if SA_strand == None:
				break
			proposed_len = delta_ref - delta_read
			if (ref_overlap < 30) and self.sv_len_pass(proposed_len, self.svlen):
				self.SA_read_supp = True
				self.SA_read_name = read.query_name
				break

	def dup_signature(self, read, mapping_quality_thr):
		# from CIGAR
		pos_list = self.cigar_dict['I']['ref_pos']
		len_list = self.cigar_dict['I']['len']
		for ind, ref_pos in enumerate(pos_list):
			if (ref_pos > self.region_buffer_left) and (ref_pos < self.region_buffer_right):
				proposed_len = len_list[ind]
				if self.sv_len_pass(proposed_len, self.svlen):
					self.CG_read_supp = True
					self.CG_read_name = read.query_name
					break

		# from supplementary alignments
		for SA_strand, delta_read, delta_ref, ref_overlap, bp1_overlap, bp2_overlap, SA_next_left, SA_next_right in self.get_sa_variables(read, mapping_quality_thr):
			if SA_strand == None:
				break
			if (self.read_strand == SA_strand):
				proposed_len = ref_overlap
				#if (float(abs(sv_len - target_svlen))/float(target_svlen) < len_ratio_tol) and \
				#	(   (float(abs(target_start - bp1_overlap))/float(target_svlen) < len_ratio_tol) or \
				#		(float(abs(target_stop - bp2_overlap))/float(target_svlen) < len_ratio_tol) ):
				if self.sv_len_pass(proposed_len, self.svlen) and \
					(   (float(abs(self.start - bp1_overlap)) < 100) or \
						(float(abs(self.stop - bp2_overlap)) < 100) ):
					self.SA_read_supp = True
					self.SA_read_name = read.query_name
					break

	def inv_signature(self, read, mapping_quality_thr):
		SA_next_left = None
		SA_next_right = None
		# we just need SA_next_left and SA_next_right after going through all the supplementary alignments
		for SA_strand, delta_read, delta_ref, ref_overlap, bp1_overlap, bp2_overlap, this_SA_next_left, this_SA_next_right in self.get_sa_variables(read, mapping_quality_thr):
			SA_next_left = this_SA_next_left
			SA_next_right = this_SA_next_right
		if SA_next_left == None:
			return

		### read has a right SA
		if (SA_next_right['SA_ref_start'] != -1):
			if self.read_strand == '+':
				breakpoint1 = self.read_ref_stop
				breakpoint2 = SA_next_right['SA_ref_stop']
			else:
				breakpoint1 = self.read_ref_start
				breakpoint2 = SA_next_right['SA_ref_start']
			min_bp_right = min(breakpoint1, breakpoint2)
			max_bp_right = max(breakpoint1, breakpoint2)
		else:
			min_bp_right = -1
			max_bp_right = 1e15
		### read has a left SA
		if (SA_next_left['SA_ref_start'] != -1):
			if self.read_strand == '+':
				breakpoint1 = self.read_ref_start
				breakpoint2 = SA_next_left['SA_ref_start']
			else:
				breakpoint1 = self.read_ref_stop
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
			proposed_len = abs(bp_2 - bp_1)
			if self.sv_len_pass(proposed_len, self.svlen):
				self.SA_read_supp = True
				self.SA_read_name = read.query_name

	def bnd_signature(self, read):
		pass

	def sv_signature(self, read, mapping_quality_thr, buffer_length):

		self.set_read_variables(read, buffer_length)

		# determine if the read is informative, meaning if any of the ends are close to the target breakends.
		### 100 is an arbitrary length
		locus_region_len = 100
		if ((self.start - locus_region_len > self.read_ref_start) and
			(self.start - locus_region_len < self.read_ref_stop)):
			self.locus_read = True
			self.locus_read_name = read.query_name
		elif ((self.start + locus_region_len > self.read_ref_start) and
			  (self.start + locus_region_len < self.read_ref_stop)):
			self.locus_read = True
			self.locus_read_name = read.query_name
		elif ((self.stop - locus_region_len > self.read_ref_start) and
			  (self.stop - locus_region_len < self.read_ref_stop)):
			self.locus_read = True
			self.locus_read_name = read.query_name
		elif ((self.stop + locus_region_len > self.read_ref_start) and
			  (self.stop + locus_region_len < self.read_ref_stop)):
			self.locus_read = True
			self.locus_read_name = read.query_name

		if not self.locus_read:
			return '', '', ''

		self.set_cigar_variables(read)

		if self.svtype == 'INS':
			self.ins_signature(read, mapping_quality_thr)
		elif self.svtype == 'DEL':
			self.del_signature(read, mapping_quality_thr)
		elif self.svtype == 'DUP':
			self.dup_signature(read, mapping_quality_thr)
		elif self.svtype == 'INV':
			self.inv_signature(read, mapping_quality_thr)
		elif self.svtype == 'BND':
			self.bnd_signature(read, mapping_quality_thr)
		else:
			raise NotImplementedError(f'this svtype is not supported for now, svtype: {self.svtype}')

		return self.locus_read_name, self.CG_read_name, self.SA_read_name

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
