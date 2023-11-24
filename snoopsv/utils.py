
class skip_class(object):

	def __init__(self, bed=None):
		self.region_list = []
		if bed:
			with open(bed, 'r') as fh:
				for line in fh.readlines():
					line = line.strip().split('\t')
					self.region_list.append({'chrom': line[0], 'start': int(line[1]), 'stop': int(line[2])})

	def skip_region(self, chrom, start, stop):

		for region in self.region_list:
			if chrom == region['chrom'] and start > region['start'] and stop < region['stop']:
				return True

		return False

	def __str__():
		ret = []
		tab = '\t'
		for region in self.region_list[:10]:
			ret.append(f"{tab.join([region['chrom'], region['start'], region['stop']])}")
		ret.append('...')
		return '\n'.join(ret)

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
	cigar_dict = {'I':{'ref_pos':[], 'len':[], 'seq_pos':[]}, 'D':{'ref_pos':[], 'seq_pos':[]}}
	cur_ref_pos = read_ref_start
	cur_seq_pos = 0
	for c in cigar_t[ind_start:ind_end+1]:
		if c[0] == 0 or c[0] == 7 or c[0] == 8: # M, =, X
			cur_ref_pos += c[1]
			cur_seq_pos += c[1]
		elif c[0] == 1: # I
			if c[1] >= ins_len_thr:
				cigar_dict['I']['ref_pos'].append(cur_ref_pos)
				cigar_dict['I']['seq_pos'].append((cur_seq_pos, cur_seq_pos+c[1]))
				cigar_dict['I']['len'].append(c[1])
			cur_seq_pos += c[1]
		elif c[0] == 2: # D
			if c[1] >= del_len_thr:
				cigar_dict['D']['ref_pos'].append((cur_ref_pos, cur_ref_pos+c[1]))
				cigar_dict['D']['seq_pos'].append(cur_seq_pos)
			cur_ref_pos += c[1]
	assert cur_ref_pos == read_ref_stop, 'something wrong with CIGAR length addition, cur_ref_pos: '+str(cur_ref_pos)+', read_ref_stop: '+str(read_ref_stop)
	assert cur_seq_pos == read_al_stop-read_al_start, 'something wrong with CIGAR length addition, cur_seq_pos: '+str(cur_seq_pos)+', read_al_stop-read_al_start: '+str(read_al_stop-read_al_start)

	return cigar_dict, ind_start, ind_end

def merge_cigar_dels(ref_pos_tuple_list, seq_pos_list, merge_fraction_thr=0.05):
	assert len(seq_pos_list) == len(ref_pos_tuple_list), f'problem with input. ref_pos_tuple_list: {ref_pos_tuple_list}, seq_pos_list: {seq_pos_list}'
	if len(ref_pos_tuple_list) < 2:
		return ref_pos_tuple_list
	ret_ref_pos_tuple_list = [x for x in ref_pos_tuple_list]
	# number of bps between two Ds are stored in the left D (in the CIGAR orientation)
	ret_nbp_within_list = [abs(x - y) for x, y in zip(seq_pos_list[:-1], seq_pos_list[1:])] + [0]
	assert len(ret_ref_pos_tuple_list) == len(ret_nbp_within_list), f'proble with ret_ref_pos_tuple_list: {ret_ref_pos_tuple_list}, ret_nbp_within_list: {ret_nbp_within_list}'
	#print(f'ref_pos_tuple_list: {ref_pos_tuple_list}')
	#print(f'seq_pos_list: {seq_pos_list}')
	#print(f'ret_ref_pos_tuple_list: {ret_ref_pos_tuple_list}')
	#print(f'ret_nbp_within_list: {ret_nbp_within_list}')
	lists_modified = True
	while lists_modified and len(ret_ref_pos_tuple_list) >= 2:
		pre_ref_pos_tuple = ret_ref_pos_tuple_list[0]
		pre_nbp_within = ret_nbp_within_list[0]
		for ind, cur_ref_pos_tuple in enumerate(ret_ref_pos_tuple_list):
			if ind == 0:
				continue
			cur_nbp_within = ret_nbp_within_list[ind]
			nbp_within = pre_nbp_within #cur_nbp_within  + pre_nbp_within
			len_merge = cur_ref_pos_tuple[1] - pre_ref_pos_tuple[0]
			if (nbp_within / len_merge) < merge_fraction_thr:
				# overwrite the left D and pop the right D
				ret_ref_pos_tuple_list[ind - 1] = (pre_ref_pos_tuple[0], cur_ref_pos_tuple[1])
				ret_nbp_within_list[ind - 1] = nbp_within + cur_nbp_within
				ret_ref_pos_tuple_list.pop(ind)
				ret_nbp_within_list.pop(ind)
				lists_modified = True
				#print(f'merging...')
				#print(f'ret_ref_pos_tuple_list: {ret_ref_pos_tuple_list}')
				#print(f'ret_nbp_within_list: {ret_nbp_within_list}')
				break
			pre_ref_pos_tuple = cur_ref_pos_tuple
			pre_nbp_within = cur_nbp_within
			lists_modified = False

	#print('Finished...')
	#print(f'ret_ref_pos_tuple_list: {ret_ref_pos_tuple_list}')
	#print(f'ret_nbp_within_list: {ret_nbp_within_list}')
	return ret_ref_pos_tuple_list

def calc_recip_overlap(s1, e1, s2, e2):
	assert (s1 <= e1) and (s2 <= e2), f'input error, s1: {s1}, e1: {e1}, s2: {s2}, e2: {e2}'
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
		if c[0] == 0 or c[0] == 7 or c[0] == 8: # Match
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
		if c[0] == 0 or c[0] == 7 or c[0] == 8: # Match
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
		if c[0] == 0 or c[0] == 7 or c[0] == 8: # Match
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
