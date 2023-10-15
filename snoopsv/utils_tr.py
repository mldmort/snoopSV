import pysam
from snoopsv.utils import get_SA_cigar_dict, get_seq_segment, get_seq_segment_supp

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

