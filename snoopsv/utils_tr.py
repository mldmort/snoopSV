import pysam
from numpy import median
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

	#Q_1 = 0
	if len(h1_list) >= r_min:
		h1_gt = str(int(median(h1_list)))
		#Q_1 =
	else:
		h1_gt = '.'

	#Q_2 = 0
	if len(h2_list) >= r_min:
		h2_gt = str(int(median(h2_list)))
		#Q_2 =
	else:
		h2_gt = '.'

	return h1_gt+'|'+h2_gt

def tr_signature_3(read, tr_start, tr_end, visited_read_set, bam, mapping_quality_thr, flanking_bp, verbose_debug):
	
	read_ref_start = read.reference_start
	read_ref_stop = read.reference_end
	read_strand = '+'
	if read.is_reverse:
		read_strand = '-'

	#if verbose_debug:
	#	print('===================================================')
	#	print('read.query_name', read.query_name)
	#	print('tr_start:', tr_start)
	#	print('tr_end:', tr_end)
	#	print('read_ref_start:', read_ref_start)
	#	print('read_ref_stop:', read_ref_stop)
	#	print('read.query_alignment_start:', read.query_alignment_start)
	#	print('read.query_alignment_end:', read.query_alignment_end)

	#### we assume pysam fetches reads from smallest ref position, so we first encounter the left of the tr region
	#### find SA alignments which span tr_end + flanking_bp. We check if the read itself spans tr_start - flanking_bp later
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
			#if verbose_debug:
			#	print('SA:', SA)
			#	print('SA_chrom:', SA_chrom)
			#	print('SA_ref_start:', SA_ref_start)
			#	print('SA_ref_stop:', SA_ref_stop)
			#	print('SA_strand:', SA_strand)
			#	print('SA_cigar:', SA_cigar)
			#	print('SA_mapq:', SA_mapq)
			if (read.reference_name == SA_chrom) and (SA_ref_start < (tr_end + flanking_bp)) and (SA_ref_stop > (tr_end + flanking_bp)) and (read_strand == SA_strand) and (SA_mapq >= mapping_quality_thr):
				has_sa_read = True
				#if verbose_debug:
				#	print('-------> found a SA read...')
				break

	sa_supp = False
	SA_reads = []
	if has_sa_read:
		key_read_self = '_'.join([str(read_ref_start), str(read_ref_stop), read.get_tag(tag="SA")])
		fh_bam = pysam.AlignmentFile(bam, 'rb')
		for sa_read in fh_bam.fetch(read.reference_name, max(0, tr_start - flanking_bp), tr_end + flanking_bp):
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
				if (read.reference_name == sa_read.reference_name) and (sa_read_ref_start < (tr_end + flanking_bp)) and (sa_read_ref_stop > (tr_end + flanking_bp)) and (read_strand == sa_read_strand) and (sa_read_mapq >= mapping_quality_thr) and (key_sa_read != key_read_self):
					SA_reads.append(sa_read)
					#if verbose_debug:
					#	print('-------> setting this read as SA read...')
					#	print('sa_read_ref_start:', sa_read_ref_start)
					#	print('sa_read_ref_stop:', sa_read_ref_stop)
					#	print('sa_read_strand:', sa_read_strand)
					#	print('sa_read_mapq:', sa_read_mapq)
					#	print('sa_read.query_alignment_start:', sa_read.query_alignment_start)
					#	print('sa_read.query_alignment_end:', sa_read.query_alignment_end)
		fh_bam.close()

	### query_alignment_start: this always starts from the left side of CIGAR, no mater + or - strand read
	### so the largest one is the right most alignment in the ref coordinates. That is the one we choose to process
	SA_reads.sort(key=lambda x: x.query_alignment_start) # sorts from small to large, the last one is the largest
	#print('SA_reads:', SA_reads)

	### now check if the read itself qualifies, and choose a mode: SA mode, CIGAR mode
	### in SA mode: read should span tr_start - flanking_bp, and there must be a qualified SA read which spans tr_end + flanking_bp
	### in CIGAR mode: read should span the whole tr region
	if (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_start - flanking_bp)) and (len(SA_reads) > 0) and (read.query_name not in visited_read_set):
		#if verbose_debug:
		#	print('in SA mode:')
		sa_supp = True
		locus_read_name = read.query_name
		sa_read = SA_reads[-1]
		read_seq_tr, blank_start, blank_stop = get_seq_segment_supp(read, sa_read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		num_bp = len(read_seq_tr) - max(0, flanking_bp - blank_start) - max(0, flanking_bp - blank_stop)
	elif (read_ref_start < (tr_start - flanking_bp)) and (read_ref_stop > (tr_end + flanking_bp)) and (read.query_name not in visited_read_set):
		#if verbose_debug:
		#	print('in CIGAR mode:')
		locus_read_name = read.query_name
		read_seq_tr, blank_start, blank_stop = get_seq_segment(read, tr_start-flanking_bp, tr_end+flanking_bp) ### tr_start and stop are 0-based
		num_bp = len(read_seq_tr) - max(0, flanking_bp - blank_start) - max(0, flanking_bp - blank_stop)
		#if verbose_debug:
		#	print(f'flanking_bp: {flanking_bp}')
		#	print(f'blank_start: {blank_start}')
		#	print(f'blank_stop: {blank_stop}')

		### align the repeat to the sequence if necessary
		#score_ind_list, raw_score_ind_list = AlignmentScore(period_seq, read_seq_tr, k_s_dict)
		#num_repeat_align = len(score_ind_list)
	else:
		#if verbose_debug:
		#	print('in Null mode:')
		locus_read_name = ''
		num_bp = -1

	#if verbose_debug:
	#	print(f'num_bp: {num_bp}')
	#	print('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
	return locus_read_name, num_bp, sa_supp
