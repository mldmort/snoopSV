import pysam
import subprocess
from snoopsv.utils import skip_region
from snoopsv.utils_tr import tr_signature_3, infer_gt_tr_phased 
from snoopsv.utils_vcf import add_header_lines

def GT_TR(tr_annot_file, vcf_in, vcf_out, contig, sample_bam_file, n_sec, i_sec, tr_span_max, r_min, verbose=1):

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
	#for line in new_header_INFO+new_header_FORMAT:
	#	header_in.add_line(line)
	add_header_lines(header_in)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	count_skip_region = 0
	for tr_annot in tr_annots_list:
		tr_chrom, tr_start, tr_end, period_len, CN, period_seq = tr_annot.split('\t')
		tr_start = int(tr_start)
		tr_end = int(tr_end)
		try:
			period_len = int(period_len)
		except:
			period_len = period_len
		try:
			CN = float(CN)
		except:
			CN = CN

		rec = fh_vcf_out.new_record(contig=tr_chrom, start=tr_start-1, stop=tr_end, alleles=('.', '.'))

		rec.info['TR_REPEAT_LEN'] = str(period_len)
		rec.info['TR_REPEAT_SEQ'] = period_seq
		rec.info['TR_REPEAT_START'] = str(tr_start)
		rec.info['TR_REPEAT_END'] = str(tr_end)
		rec.info['TR_REPEAT_CN'] = str(CN)
		rec.info['TR_ANNOT'] = True

		if skip_region(tr_chrom, tr_start, tr_end):
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
				#if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr) and (read.query_name == '257dfc21-24cf-4751-891a-107cc685e9ad'):
					#print('working on this:', read.query_name)
					locus_read_name, num_repeat_al, num_repeat_ln, num_bp = tr_signature_3(read, tr_start, tr_end, period_len, CN, period_seq, k_s_dict, visited_read_set, bam_file, mapping_quality_thr)
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
							tr_supp_ln['H'+str(HP)].append(num_repeat_ln)
						else:
							tr_supp_ln['H0'].append(num_repeat_ln)
					if num_bp >= 0:
						if read.has_tag('HP'):
							HP = read.get_tag(tag='HP')
							tr_supp_bp['H'+str(HP)].append(num_bp)
						else:
							tr_supp_bp['H0'].append(num_bp)
			fh_bam.close()

			#tr_GT_ln = infer_gt_tr_phased(tr_supp_ln, r_min)
			tr_GT_bp = infer_gt_tr_phased(tr_supp_bp, r_min)
			#rec.samples[sample]['GT_TR_LN'] = tr_GT_ln
			rec.samples[sample]['GT_TR_BP'] = tr_GT_bp

			tr_len = tr_end - tr_start
			bp_list = tr_GT_bp.split('|')
			bp_dv_list = []
			if bp_list[0]=='.':
				bp_dv_list.append('.')
			else:
				bp_dv_list.append(str( int(bp_list[0]) - tr_len ))
			if bp_list[1]=='.':
				bp_dv_list.append('.')
			else:
				bp_dv_list.append(str( int(bp_list[1]) - tr_len ))
			rec.samples[sample]['BP_DV'] = '|'.join(bp_dv_list)

			#temp = '|'.join([str(x) for x in tr_supp_ln['H1']])
			#if temp == '':
			#	temp = '.'
			#rec.samples[sample]['CN_TR_LN_H1'] = temp

			#temp = '|'.join([str(x) for x in tr_supp_ln['H2']])
			#if temp == '':
			#	temp = '.'
			#rec.samples[sample]['CN_TR_LN_H2'] = temp

			#temp = '|'.join([str(x) for x in tr_supp_ln['H0']])
			#if temp == '':
			#	temp = '.'
			#rec.samples[sample]['CN_TR_LN_H0'] = temp

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

