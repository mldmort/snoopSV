from math import log10
from snoopsv.utils import get_cigar_dict, calc_recip_overlap, merge_cigar_dels

class sv_class:

	def __init__(self, rec, len_ratio_tol, ins_len_thr, del_len_thr, del_recip_overlap_thr,
				 ins_recip_overlap_thr, dup_recip_overlap_thr, inv_recip_overlap_thr, bnd_pos_tol):

		self.len_ratio_tol = len_ratio_tol
		self.ins_len_thr = ins_len_thr
		self.del_len_thr = del_len_thr
		self.del_recip_overlap_thr = del_recip_overlap_thr
		self.ins_recip_overlap_thr = ins_recip_overlap_thr
		self.dup_recip_overlap_thr = dup_recip_overlap_thr
		self.inv_recip_overlap_thr = inv_recip_overlap_thr
		self.bnd_pos_tol = bnd_pos_tol

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
			except KeyError:
				raise KeyError(f'SVLEN is not defined for id: {rec.id}. Check the VCF file!')
			except TypeError:
				if isinstance(rec.info['SVLEN'], tuple):
					assert len(rec.info['SVLEN']) == 1
					self.svlen = rec.info['SVLEN'][0]
				else:
					raise TypeError(f'problem with SVLEN. rec.info["SVLEN"]: {rec.info["SVLEN"]}, rec.id: {rec.id}')
		else:
			self.svlen = None

		self.chr2 = self.chrom # for non BND svtypes
		if self.svtype == 'BND':
			if 'CHR2' in rec.info:
				self.chr2 = rec.info['CHR2']
			if '[' in rec.alts[0]:
				chr2_pos2 = rec.alts[0].split('[')[1] # always second member, both N[... and [...
				chr2_pos2 = chr2_pos2.split(':')
				self.chr2 = chr2_pos2[0]
				self.pos2 = int(chr2_pos2[1])
			elif ']' in rec.alts[0]:
				chr2_pos2 = rec.alts[0].split(']')[1] # always second member, both N]... and ]...
				chr2_pos2 = chr2_pos2.split(':')
				self.chr2 = chr2_pos2[0]
				self.pos2 = int(chr2_pos2[1])
			else:
				raise NameError(f'[ or ] not found in alt column of BND svtype, for sv id: {self.id}')
			#assert self.chr2 != self.chrom, f'problem with parsing chr2: {self.chr2}, chrom: {self.chrom}, sv id: {self.id}'
		else:
			self.pos2 = None

	def __str__(self):
		ret = (f'++++++++++++++++++++\n'
			   f'id: {self.id}\n'
			   f'chrom: {self.chrom}\n'
			   f'start: {self.start}\n'
			   f'stop: {self.stop}\n'
			   f'svtype: {self.svtype}\n'
			   f'svlen: {self.svlen}\n'
			   f'chr2: {self.chr2}\n'
			   f'pos2: {self.pos2}\n'
			   f'++++++++++++++++++++')
		return ret

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

	# this gives the left and right supplementary alignments if present
	def get_sa_variables(self, read, mapping_quality_thr):
		SA_next_right = {'SA_chrom': None, 'SA_strand': None, 'SA_ref_start': -1, 'SA_ref_stop': -1, 'SA_read_start': 1e15, 'SA_read_stop': 1e15}
		SA_next_left = {'SA_chrom': None, 'SA_strand': None, 'SA_ref_start': -1, 'SA_ref_stop': -1, 'SA_read_start': -1, 'SA_read_stop': -1}
		if read.has_tag('SA'):
			SA_tag = read.get_tag(tag='SA')
			SA_list = SA_tag.split(';')[:-1]
			for SA in SA_list:
				SA = SA.split(',')
				SA_chrom = SA[0]
				SA_ref_start = int(SA[1])
				SA_strand = SA[2]
				SA_cigar = SA[3]
				SA_mapq = float(SA[4])

				# don't consider this SA if:
				if (SA_mapq < mapping_quality_thr):
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

				delta_read = None
				delta_ref = None
				ref_overlap = None
				bp1_delta_ref = None
				bp2_delta_ref = None
				bp1_overlap = None
				bp2_overlap = None
				bp1_inv = None
				bp2_inv = None
				bp_delta_read = None
				if (self.read_al_start < SA_read_start):
					delta_read = SA_read_start - self.read_al_stop
					if (self.read_strand == '+'):
						bp_delta_read = self.read_ref_stop
						# +++read+++> <---SA--- +++>  or  <---SA--- +++read+++> <---
						bp1_inv = min(self.read_ref_stop, SA_ref_stop)
						bp2_inv = max(self.read_ref_stop, SA_ref_stop)
					else:
						bp_delta_read = self.read_ref_start
						# <--- +++SA+++> <---read---  or  +++SA+++> <---read--- +++>
						bp1_inv = min(self.read_ref_start, SA_ref_start)
						bp2_inv = max(self.read_ref_start, SA_ref_start)
					if (self.read_strand == '+') and (SA_strand == '+'):
						delta_ref = SA_ref_start - self.read_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = SA_ref_start
						bp2_overlap = self.read_ref_stop
						bp1_delta_ref = self.read_ref_stop
						bp2_delta_ref = SA_ref_start
					elif (self.read_strand == '-') and (SA_strand == '-'):
						delta_ref = self.read_ref_start - SA_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = self.read_ref_start
						bp2_overlap = SA_ref_stop
						bp1_delta_ref = SA_ref_stop
						bp2_delta_ref = self.read_ref_start
				else:
					delta_read = self.read_al_start - SA_read_stop
					if (self.read_strand == '+'):
						bp_delta_read = self.read_ref_start
						# <--- +++read+++> <---SA---  or  +++> <---SA--- +++read+++>
						bp1_inv = min(self.read_ref_start, SA_ref_start)
						bp2_inv = max(self.read_ref_start, SA_ref_start)
					else:
						bp_delta_read = self.read_ref_stop
						# +++SA+++> <---read--- +++>  or  <---read--- +++SA+++> <---
						bp1_inv = min(self.read_ref_stop, SA_ref_stop)
						bp2_inv = max(self.read_ref_stop, SA_ref_stop)
					if (self.read_strand == '+') and (SA_strand == '+'):
						delta_ref = self.read_ref_start - SA_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = self.read_ref_start
						bp2_overlap = SA_ref_stop
						bp1_delta_ref = SA_ref_stop
						bp2_delta_ref = self.read_ref_start
					elif (self.read_strand == '-') and (SA_strand == '-'):
						delta_ref = SA_ref_start - self.read_ref_stop
						ref_overlap = -1 * delta_ref
						bp1_overlap = SA_ref_start
						bp2_overlap = self.read_ref_stop
						bp1_delta_ref = self.read_ref_stop
						bp2_delta_ref = SA_ref_start

				if (SA_read_start > self.read_al_start) and (SA_read_start < SA_next_right['SA_read_start']): # SA is right of read
					SA_next_right['SA_read_start'] = SA_read_start
					SA_next_right['SA_read_stop'] = SA_read_stop
					SA_next_right['SA_ref_start'] = SA_ref_start
					SA_next_right['SA_ref_stop'] = SA_ref_stop
					SA_next_right['SA_chrom'] = SA_chrom
					SA_next_right['SA_strand'] = SA_strand
					SA_next_right['delta_read'] = delta_read
					SA_next_right['delta_ref'] = delta_ref
					SA_next_right['ref_overlap'] = ref_overlap
					SA_next_right['bp1_overlap'] = bp1_overlap
					SA_next_right['bp2_overlap'] = bp2_overlap
					SA_next_right['bp1_delta_ref'] = bp1_delta_ref
					SA_next_right['bp2_delta_ref'] = bp2_delta_ref
					SA_next_right['bp1_inv'] = bp1_inv
					SA_next_right['bp2_inv'] = bp2_inv
					SA_next_right['bp_delta_read'] = bp_delta_read
				if (SA_read_start < self.read_al_start) and (SA_read_start > SA_next_left['SA_read_start']): # SA is left of read
					SA_next_left['SA_read_start'] = SA_read_start
					SA_next_left['SA_read_stop'] = SA_read_stop
					SA_next_left['SA_ref_start'] = SA_ref_start
					SA_next_left['SA_ref_stop'] = SA_ref_stop
					SA_next_left['SA_chrom'] = SA_chrom
					SA_next_left['SA_strand'] = SA_strand
					SA_next_left['delta_read'] = delta_read
					SA_next_left['delta_ref'] = delta_ref
					SA_next_left['ref_overlap'] = ref_overlap
					SA_next_left['bp1_overlap'] = bp1_overlap
					SA_next_left['bp2_overlap'] = bp2_overlap
					SA_next_left['bp1_delta_ref'] = bp1_delta_ref
					SA_next_left['bp2_delta_ref'] = bp2_delta_ref
					SA_next_left['bp1_inv'] = bp1_inv
					SA_next_left['bp2_inv'] = bp2_inv
					SA_next_left['bp_delta_read'] = bp_delta_read

				#sa_params = {'SA_strand': SA_strand,
				#			'SA_chrom': SA_chrom,
				#			'delta_read': delta_read,
				#			'delta_ref': delta_ref,
				#			'ref_overlap': ref_overlap,
				#			'bp1_overlap': bp1_overlap,
				#			'bp2_overlap': bp2_overlap,
				#			'SA_next_left': SA_next_left,
				#			'SA_next_right': SA_next_right,
				#			'SA_ref_start': SA_ref_start,
				#			'SA_ref_stop':  SA_ref_stop,
				#			'SA_read_start': SA_read_start,
				#			'SA_read_stop': SA_read_stop,
				#			'bp1_delta_ref': bp1_delta_ref,
				#			'bp2_delta_ref': bp2_delta_ref,
				#			}
				#yield sa_params
		return SA_next_left, SA_next_right

	def ins_signature(self, read, mapping_quality_thr):
		'''
		this method is searching for insertion signature and
		sets these parameters:
		- self.CG_read_supp
		- self.CG_read_name
		- self.SA_read_supp
		- self.SA_read_name
		'''
		#print('===================================')
		#print(f'working on read: {read.query_name}')
		#print('SV:')
		#print(self)
		#print(f'self.read_ref_start: {self.read_ref_start}')
		#print(f'self.read_ref_stop: {self.read_ref_stop}')
		#print(f'self.read_al_start: {self.read_al_start}')
		#print(f'self.read_al_stop: {self.read_al_stop}')
		#print(f'self.cigar_dict: {self.cigar_dict}')
		# from CIGAR
		pos_list = self.cigar_dict['I']['ref_pos']
		len_list = self.cigar_dict['I']['len']
		for ind, ref_pos in enumerate(pos_list):
			proposed_len = len_list[ind]
			bp1 = ref_pos - int(proposed_len/2.)
			bp2 = ref_pos + int(proposed_len/2.)
			target_bp1 = self.start - int(self.svlen/2.)
			target_bp2 = self.start + int(self.svlen/2.)
			recip_overlap = calc_recip_overlap(bp1, bp2, target_bp1, target_bp2, self)
			#print(f'CIGAR proposed_len: {proposed_len}')
			#print(f'bp1: {bp1}')
			#print(f'bp2: {bp2}')
			#print(f'target_bp1: {target_bp1}')
			#print(f'target_bp2: {target_bp2}')
			#print(f'recip_overlap: {recip_overlap}')
			if recip_overlap >= self.ins_recip_overlap_thr and self.sv_len_pass(proposed_len, self.svlen):
				self.CG_read_supp = True
				self.CG_read_name = read.query_name
				break

		#print(f'self.CG_read_supp: {self.CG_read_supp}')

		SA_next_left, SA_next_right = self.get_sa_variables(read, mapping_quality_thr)

		#print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
		#print(f'SA_next_left: {SA_next_left}')
		# if we have an SA on the left with the same strand
		if (SA_next_left['SA_chrom'] == self.read_chrom) and (SA_next_left['SA_strand'] == self.read_strand):
			delta_read = SA_next_left['delta_read']
			bp_delta_read = SA_next_left['bp_delta_read']
			if delta_read > 0: # just to make sure
				bp1 = bp_delta_read - int(delta_read/2.)
				bp2 = bp_delta_read + int(delta_read/2.)
				target_bp1 = self.start - int(self.svlen/2.)
				target_bp2 = self.start + int(self.svlen/2.)
				recip_overlap = calc_recip_overlap(bp1, bp2, target_bp1, target_bp2, self)
				#print(f'bp1: {bp1}')
				#print(f'bp2: {bp2}')
				#print(f'target_bp1: {target_bp1}')
				#print(f'target_bp2: {target_bp2}')
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.ins_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

			# check if the read has a DUP signal
			ref_overlap = SA_next_left['ref_overlap']
			if ref_overlap > 0:
				bp1 = bp_delta_read - int(ref_overlap/2.)
				bp2 = bp_delta_read + int(ref_overlap/2.)
				target_bp1 = self.start - int(self.svlen/2.)
				target_bp2 = self.start + int(self.svlen/2.)
				recip_overlap = calc_recip_overlap(bp1, bp2, target_bp1, target_bp2, self)
				#print(f'DUP check: ref_overlap: {ref_overlap}')
				#print(f'DUP check: bp1: {bp1}')
				#print(f'DUP check: bp2: {bp2}')
				#print(f'DUP check: target_bp1: {target_bp1}')
				#print(f'DUP check: target_bp2: {target_bp2}')
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.ins_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		#print(f'SA_next_right: {SA_next_right}')
		# if we have an SA on the right with the same strand
		if (SA_next_right['SA_chrom'] == self.read_chrom) and (SA_next_right['SA_strand'] == self.read_strand):
			delta_read = SA_next_right['delta_read']
			bp_delta_read = SA_next_right['bp_delta_read']
			if delta_read > 0: # just to make sure
				bp1 = bp_delta_read - int(delta_read/2.)
				bp2 = bp_delta_read + int(delta_read/2.)
				target_bp1 = self.start - int(self.svlen/2.)
				target_bp2 = self.start + int(self.svlen/2.)
				recip_overlap = calc_recip_overlap(bp1, bp2, target_bp1, target_bp2, self)
				#print(f'bp1: {bp1}')
				#print(f'bp2: {bp2}')
				#print(f'target_bp1: {target_bp1}')
				#print(f'target_bp2: {target_bp2}')
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.ins_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

			# check if the read has a DUP signal
			ref_overlap = SA_next_right['ref_overlap']
			if ref_overlap > 0:
				bp1 = bp_delta_read - int(ref_overlap/2.)
				bp2 = bp_delta_read + int(ref_overlap/2.)
				target_bp1 = self.start - int(self.svlen/2.)
				target_bp2 = self.start + int(self.svlen/2.)
				recip_overlap = calc_recip_overlap(bp1, bp2, target_bp1, target_bp2, self)
				#print(f'DUP check: ref_overlap: {ref_overlap}')
				#print(f'DUP check: bp1: {bp1}')
				#print(f'DUP check: bp2: {bp2}')
				#print(f'DUP check: target_bp1: {target_bp1}')
				#print(f'DUP check: target_bp2: {target_bp2}')
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.ins_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name
		#print(f'self.SA_read_supp: {self.SA_read_supp}')

	def del_signature(self, read, mapping_quality_thr):
		'''
		this method is searching for deletion signature and
		sets these parameters:
		- self.CG_read_supp
		- self.CG_read_name
		- self.SA_read_supp
		- self.SA_read_name
		'''
		#print('===================================')
		#print(f'working on read: {read.query_name}')
		#print('SV:')
		#print(self)
		#print(f'self.read_ref_start: {self.read_ref_start}')
		#print(f'self.read_ref_stop: {self.read_ref_stop}')
		#print(f'self.read_al_start: {self.read_al_start}')
		#print(f'self.read_al_stop: {self.read_al_stop}')
		# from CIGAR
		ref_pos_list = merge_cigar_dels(self.cigar_dict['D']['ref_pos'], self.cigar_dict['D']['seq_pos'])
		#print(f'self.cigar_dict: {self.cigar_dict}')
		#print(f'ref_pos_list: {ref_pos_list}')
		for ref_pos_t in ref_pos_list:
			ref_pos_start = ref_pos_t[0]
			ref_pos_stop = ref_pos_t[1]
			#print(f'ref_pos_start: {ref_pos_start}')
			#print(f'ref_pos_stop: {ref_pos_stop}')
			if (ref_pos_start > self.region_buffer_left) and (ref_pos_stop < self.region_buffer_right):
				recip_overlap = calc_recip_overlap(ref_pos_start, ref_pos_stop, self.start, self.stop, self)
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.del_recip_overlap_thr:
					self.CG_read_supp = True
					self.CG_read_name = read.query_name
					break

		SA_next_left, SA_next_right = self.get_sa_variables(read, mapping_quality_thr)

		#print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
		#print(f'SA_next_left: {SA_next_left}')
		# if we have an SA on the left with the same strand
		if (SA_next_left['SA_chrom'] == self.read_chrom) and (SA_next_left['SA_strand'] == self.read_strand):
			bp1_delta_ref = SA_next_left['bp1_delta_ref']
			bp2_delta_ref = SA_next_left['bp2_delta_ref']
			if bp1_delta_ref < bp2_delta_ref:
				recip_overlap = calc_recip_overlap(bp1_delta_ref, bp2_delta_ref, self.start, self.stop, self)
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.del_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		#print(f'SA_next_right: {SA_next_right}')
		# if we have an SA on the right with the same strand
		if (SA_next_right['SA_chrom'] == self.read_chrom) and (SA_next_right['SA_strand'] == self.read_strand):
			bp1_delta_ref = SA_next_right['bp1_delta_ref']
			bp2_delta_ref = SA_next_right['bp2_delta_ref']
			if bp1_delta_ref < bp2_delta_ref:
				recip_overlap = calc_recip_overlap(bp1_delta_ref, bp2_delta_ref, self.start, self.stop, self)
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.del_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		#print(f'self.CG_read_supp: {self.CG_read_supp}')
		#print(f'self.SA_read_supp: {self.SA_read_supp}')

	def dup_signature(self, read, mapping_quality_thr):
		'''
		this method is searching for duplication signature and
		sets these parameters:
		- self.CG_read_supp
		- self.CG_read_name
		- self.SA_read_supp
		- self.SA_read_name
		'''
		#print('===================================')
		#print(f'working on read: {read.query_name}')
		#print('SV:')
		#print(self)
		#print(f'self.read_ref_start: {self.read_ref_start}')
		#print(f'self.read_ref_stop: {self.read_ref_stop}')
		#print(f'self.read_al_start: {self.read_al_start}')
		#print(f'self.read_al_stop: {self.read_al_stop}')
		#print(f'self.cigar_dict: {self.cigar_dict}')
		# from CIGAR, signal from insertion calls
		pos_list = self.cigar_dict['I']['ref_pos']
		len_list = self.cigar_dict['I']['len']
		for ind, ref_pos in enumerate(pos_list):
			proposed_len = len_list[ind]
			bp1_op1 = ref_pos - proposed_len
			bp2_op1 = ref_pos
			bp1_op2 = ref_pos
			bp2_op2 = ref_pos + proposed_len
			recip_overlap_op1 = calc_recip_overlap(bp1_op1, bp2_op1, self.start, self.stop, self)
			recip_overlap_op2 = calc_recip_overlap(bp1_op2, bp2_op2, self.start, self.stop, self)
			#print(f'CIGAR proposed_len: {proposed_len}')
			#print(f'bp1_op1: {bp1_op1}')
			#print(f'bp2_op1: {bp2_op1}')
			#print(f'recip_overlap_op1: {recip_overlap_op1}')
			#print(f'bp1_op2: {bp1_op2}')
			#print(f'bp2_op2: {bp2_op2}')
			#print(f'recip_overlap_op2: {recip_overlap_op2}')
			if recip_overlap_op1 >= self.dup_recip_overlap_thr or recip_overlap_op2 >= self.dup_recip_overlap_thr:
				self.CG_read_supp = True
				self.CG_read_name = read.query_name
				break

		#print(f'self.CG_read_supp: {self.CG_read_supp}')

		SA_next_left, SA_next_right = self.get_sa_variables(read, mapping_quality_thr)

		#print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
		#print(f'SA_next_left: {SA_next_left}')
		# if we have an SA on the left with the same strand
		if (SA_next_left['SA_chrom'] == self.read_chrom) and (SA_next_left['SA_strand'] == self.read_strand):
			ref_overlap = SA_next_left['ref_overlap']
			bp1_overlap = SA_next_left['bp1_overlap']
			bp2_overlap = SA_next_left['bp2_overlap']
			#print(f'bp1_overlap: {bp1_overlap}')
			#print(f'bp2_overlap: {bp2_overlap}')
			#print(f'ref_overlap: {ref_overlap}')
			if ref_overlap > 0:
				recip_overlap = calc_recip_overlap(bp1_overlap, bp2_overlap, self.start, self.stop, self)
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.dup_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		#print(f'SA_next_right: {SA_next_right}')
		# if we have an SA on the right with the same strand
		if (SA_next_right['SA_chrom'] == self.read_chrom) and (SA_next_right['SA_strand'] == self.read_strand):
			ref_overlap = SA_next_right['ref_overlap']
			bp1_overlap = SA_next_right['bp1_overlap']
			bp2_overlap = SA_next_right['bp2_overlap']
			#print(f'bp1_overlap: {bp1_overlap}')
			#print(f'bp2_overlap: {bp2_overlap}')
			#print(f'ref_overlap: {ref_overlap}')
			if ref_overlap > 0:
				recip_overlap = calc_recip_overlap(bp1_overlap, bp2_overlap, self.start, self.stop, self)
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.dup_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		#print(f'self.SA_read_supp: {self.SA_read_supp}')

	def inv_signature(self, read, mapping_quality_thr):
		'''
		this method is searching for inversion signature and
		sets these parameters:
		- self.CG_read_supp
		- self.CG_read_name
		- self.SA_read_supp
		- self.SA_read_name
		'''
		#print('===================================')
		#print(f'working on read: {read.query_name}')
		#print('SV:')
		#print(self)
		#print(f'self.read_ref_start: {self.read_ref_start}')
		#print(f'self.read_ref_stop: {self.read_ref_stop}')
		#print(f'self.read_al_start: {self.read_al_start}')
		#print(f'self.read_al_stop: {self.read_al_stop}')

		SA_next_left, SA_next_right = self.get_sa_variables(read, mapping_quality_thr)

		has_right_SA = False
		has_left_SA = False

		#print('xxxxxxxxxxxxxxxxxxxxxxxxxx')
		#print(f'SA_next_right: {SA_next_right}')
		# read has a right SA
		if (SA_next_right['SA_chrom'] == self.read_chrom) and (SA_next_right['SA_strand'] != self.read_strand):
			has_right_SA = True
			bp1_inv_right = SA_next_right['bp1_inv']
			bp2_inv_right = SA_next_right['bp2_inv']
			#print(f'bp1_inv_right: {bp1_inv_right}')
			#print(f'bp2_inv_right: {bp2_inv_right}')

		#print(f'SA_next_left: {SA_next_left}')
		# read has a left SA
		if (SA_next_left['SA_chrom'] == self.read_chrom) and (SA_next_left['SA_strand'] != self.read_strand):
			has_left_SA = True
			bp1_inv_left = SA_next_left['bp1_inv']
			bp2_inv_left = SA_next_left['bp2_inv']
			#print(f'bp1_inv_left: {bp1_inv_left}')
			#print(f'bp2_inv_left: {bp2_inv_left}')

		if has_right_SA and has_left_SA:
			bp1 = max(bp1_inv_right, bp1_inv_left)
			bp2 = min(bp2_inv_right, bp2_inv_left)
			if bp1 < bp2:
				recip_overlap = calc_recip_overlap(bp1, bp2, self.start, self.stop, self)
				#print(f'bp1: {bp1}')
				#print(f'bp2: {bp2}')
				#print(f'recip_overlap: {recip_overlap}')
				if recip_overlap >= self.inv_recip_overlap_thr:
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		if has_right_SA:
			bp1 = bp1_inv_right
			bp2 = bp2_inv_right
			recip_overlap = calc_recip_overlap(bp1, bp2, self.start, self.stop, self)
			#print(f'right bp1: {bp1}')
			#print(f'right bp2: {bp2}')
			#print(f'right recip_overlap: {recip_overlap}')
			if recip_overlap >= self.inv_recip_overlap_thr:
				self.SA_read_supp = True
				self.SA_read_name = read.query_name

		if has_left_SA:
			bp1 = bp1_inv_left
			bp2 = bp2_inv_left
			recip_overlap = calc_recip_overlap(bp1, bp2, self.start, self.stop, self)
			#print(f'left bp1: {bp1}')
			#print(f'left bp2: {bp2}')
			#print(f'left recip_overlap: {recip_overlap}')
			if recip_overlap >= self.inv_recip_overlap_thr:
				self.SA_read_supp = True
				self.SA_read_name = read.query_name

		#print(f'self.SA_read_supp: {self.SA_read_supp}')

	def bnd_signature(self, read, mapping_quality_thr):
		'''
		this method is searching for breakend signature and
		sets these parameters:
		- self.CG_read_supp
		- self.CG_read_name
		- self.SA_read_supp
		- self.SA_read_name
		'''
		SA_next_left, SA_next_right = self.get_sa_variables(read, mapping_quality_thr)

		# read has a right SA with right mapped chrom
		if (SA_next_right['SA_ref_start'] != -1 and SA_next_right['SA_chrom'] == self.chr2):
			if self.read_strand == '+':
				breakpoint1 = self.read_ref_stop
			else:
				breakpoint1 = self.read_ref_start
			if SA_next_right['SA_strand'] == '+':
				breakpoint2 = SA_next_right['SA_ref_start']
			else:
				breakpoint2 = SA_next_right['SA_ref_stop']
			if (abs(breakpoint1 - self.start) < self.bnd_pos_tol and
				abs(breakpoint2 - self.pos2) < self.bnd_pos_tol):
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

		# read has a left SA with right mapped chrom
		if (SA_next_left['SA_ref_start'] != -1  and SA_next_left['SA_chrom'] == self.chr2):
			if self.read_strand == '+':
				breakpoint1 = self.read_ref_start
			else:
				breakpoint1 = self.read_ref_stop
			if SA_next_left['SA_strand'] == '+':
				breakpoint2 = SA_next_left['SA_ref_stop']
			else:
				breakpoint2 = SA_next_left['SA_ref_start']
			if (abs(breakpoint1 - self.start) < self.bnd_pos_tol and
				abs(breakpoint2 - self.pos2) < self.bnd_pos_tol):
					self.SA_read_supp = True
					self.SA_read_name = read.query_name

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
	
	if (DR == 0) and (DV == 0):
		return './.', 0, 0, 0, 0, 0
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
