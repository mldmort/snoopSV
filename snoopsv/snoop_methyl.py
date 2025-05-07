import pysam
import subprocess
from snoopsv.utils import skip_class
from snoopsv.utils_methyl import methylation_signature
from snoopsv.utils_vcf import add_header_lines_methyl
import sys

def process_methylation(annot_file, vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, mapping_quality_thr, buffer_length, region_span_max, r_min, annot_columns, include_columns, exclude_columns, header_lines, info_prefix, skip_bed, flanking_bp, verbose=1):

	verbose_debug = True

	annot_has_header = False
	annot_columns_list = []
	extra_fields_list = []
	extra_fields_idx_list = []
	annot_len_dev_column_idx = None

	# check if the annotation file has a header
	out = subprocess.run(['head', '-n', '1', annot_file], check=True, capture_output=True, text=True).stdout
	temp_first_line = [x.strip() for x in out.split('\t')] # need to strip to get rid of \n for the last one
	if temp_first_line[0].startswith('#'):
		annot_has_header = True
		annot_columns_list = list(temp_first_line)

	# if annotation columns are provided use them instead of the annotation file header
	if annot_columns:
		annot_columns_list = annot_columns.split(',')
		if annot_has_header:
			assert len(temp_first_line) == len(annot_columns_list), f'length of annotation header does not match the provided columns. annotation file: {temp_first_line}, columns: {annot_columns_list}'

	# include/exclude annotation columns
	if include_columns != 'all':
		include_columns_list = include_columns.split(',')
	else:
		include_columns_list = list(annot_columns_list[3:])

	exclude_columns_list = []
	if exclude_columns:
		exclude_columns_list = exclude_columns.split(',')

	# select columns to add to INFO
	for idx, col in enumerate(annot_columns_list):
		if col in include_columns_list and col not in exclude_columns_list:
			extra_fields_list.append(col)
			extra_fields_idx_list.append(idx)

	print(f'annotation file columns: {", ".join(annot_columns_list)}')
	print(f'annotation columns added to the INFO: {", ".join(extra_fields_list)}')

	if contig:
		command = ['awk', 'BEGIN{FS="\t";OFS="\t"}$1=="'+contig+'"', annot_file]
	else:
		command = ['awk', 'BEGIN{FS="\t";OFS="\t"}$1!~"#"', annot_file]
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
	annot_list = out.split('\n')[:-1] # the last one is always an empty string
	#print('annot_list[0]:', annot_list[0])
	#print('annot_list[-1]:', annot_list[-1])

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	extra_fields_toheader_list = extra_fields_list
	add_header_lines_methyl(header_in, info_prefix, extra_fields_toheader_list)
	if header_lines:
		with open(header_lines, 'r') as fh:
			for line in fh.readlines():
				header_in.add_line(line.strip())

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)
	if sample not in fh_vcf_out.header.samples:
		print(f'Adding sample: {sample} to the output VCF')
		fh_vcf_out.header.add_sample(sample)

	fh_bam = pysam.AlignmentFile(bam, 'rb')

	count_skip_region = 0
	sys.stdout.flush()
	for annot in annot_list:
		annot_split = annot.split('\t')
		chrom = annot_split[0]
		start = int(annot_split[1])
		end = int(annot_split[2])
		annot_id = '_'.join([info_prefix, chrom, str(start), str(end)])

		rec = fh_vcf_out.new_record(contig=chrom, start=start, stop=end, alleles=('.', '.'), id=annot_id)

		for extra_field_idx, extra_field in zip(extra_fields_idx_list, extra_fields_list):
			rec.info[info_prefix + '_' + extra_field] = annot_split[extra_field_idx]
		#rec.info[info_prefix + '_' + 'ANNOT'] = True
		rec.stop = end

		skip = skip_class(skip_bed)
		if skip.skip_region(chrom, start, end):
			count_skip_region += 1
			rec.info['SKIP_REGION'] = True
			fh_vcf_out.write(rec)
			continue

		visited_read_set = set()
		methyl_prob_dict = {'H1':[], 'H2':[], 'H0':[]}
		num_bp_dict = {'H1':[], 'H2':[], 'H0':[]}
		read_name_dict = {'H1':[], 'H2':[], 'H0':[]}
		phase_set_dict = {'H1':[], 'H2':[]}

		#sa_supp_any = False
		for read in fh_bam.fetch(chrom, max(0, start - buffer_length), end + buffer_length):
			if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
			### ont sample
			#if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr) and (read.query_name == '3014939c-d797-4fa6-b669-8eec48317689'):
				#print('working on this:', read.query_name)
				methyl_probs, num_bp = methylation_signature(read, start, end, flanking_bp, verbose_debug)
				if len(methyl_probs) > 0:
					methyl_probs_str = '|'.join([str(x) for x in methyl_probs])
					if read.has_tag('HP'):
						HP = read.get_tag(tag='HP')
						methyl_prob_dict['H'+str(HP)].append(methyl_probs_str)
						num_bp_dict['H'+str(HP)].append(num_bp)
						read_name_dict['H'+str(HP)].append(read.query_name)
						if read.has_tag('PS'):
							PS = read.get_tag(tag='PS')
						else:
							PS = '.'
						phase_set_dict['H'+str(HP)].append(str(PS))
					else:
						methyl_prob_dict['H0'].append(methyl_probs_str)
						num_bp_dict['H0'].append(num_bp)
						read_name_dict['H0'].append(read.query_name)

		length = end - start

		temp = ','.join([str(x) for x in methyl_prob_dict['H1']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['METHYL_H1'] = temp

		temp = ','.join([str(x) for x in methyl_prob_dict['H2']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['METHYL_H2'] = temp

		temp = ','.join([str(x) for x in methyl_prob_dict['H0']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['METHYL_H0'] = temp

		temp = ','.join([str(x - length) if x!='.' else '.' for x in num_bp_dict['H1']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['BP_DEV_H1'] = temp

		temp = ','.join([str(x - length) if x!='.' else '.' for x in num_bp_dict['H2']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['BP_DEV_H2'] = temp

		temp = ','.join([str(x - length) if x!='.' else '.' for x in num_bp_dict['H0']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['BP_DEV_H0'] = temp

		temp = ','.join([x for x in read_name_dict['H1']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['RN_H1'] = temp

		temp = ','.join([x for x in read_name_dict['H2']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['RN_H2'] = temp

		temp = ','.join([x for x in read_name_dict['H0']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['RN_H0'] = temp

		rec.samples[sample]['N_H1'] = len(num_bp_dict['H1'])
		rec.samples[sample]['N_H2'] = len(num_bp_dict['H2'])
		rec.samples[sample]['N_H0'] = len(num_bp_dict['H0'])

		temp = ','.join([x for x in phase_set_dict['H1']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['PS_H1'] = temp

		temp = ','.join([x for x in phase_set_dict['H2']])
		if temp == '':
			temp = '.'
		rec.samples[sample]['PS_H2'] = temp

		fh_vcf_out.write(rec)

	if verbose == 1:
		print('Finished all variants')
		print('count_skip_region:', count_skip_region)

	fh_bam.close()
	fh_vcf_in.close()
	fh_vcf_out.close()
