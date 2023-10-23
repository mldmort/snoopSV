import pysam
import subprocess
from snoopsv.utils import skip_class
from snoopsv.utils_sv import sv_class, infer_gt_sv, get_phased_gt
from snoopsv.utils_vcf import add_header_lines
from pathlib import Path
import itertools
import sys

def GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, skip_bed, mapping_quality_thr, buffer_length, p_err, len_ratio_tol, ins_len_thr, del_len_thr, del_recip_overlap_thr, bnd_pos_tol, verbose=1, include_svtype=None, exclude_svtype=None):

	# count the number of variants to be processed
	if contig:
		command = ('bcftools query -r ' + contig + ' -f %CHROM\\n ' + vcf_in).split()
	else:
		command = ('bcftools query -f %CHROM\\n ' + vcf_in).split()
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
	add_header_lines(fh_vcf_in.header)
	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=fh_vcf_in.header)
	# after the next few lines output header have sample but input header does not
	if sample not in fh_vcf_out.header.samples:
		print(f'Adding sample: {sample} to the output VCF')
		fh_vcf_out.header.add_sample(sample)

	fh_bam = pysam.AlignmentFile(bam, 'rb')

	skip = skip_class(skip_bed)

	count_skip_region = 0
	count_skip_sec = 0
	count_skip_svtype = 0
	sys.stdout.flush()
	for i_rec, rec in enumerate(fh_vcf_in.fetch(contig=contig)):
		if (i_rec < i_rec_start) or (i_rec >= i_rec_end):
			count_skip_sec += 1
			continue
		# if sample is not in the input header we need to create a new record to add the sample.
		# we cannot add a new sample to a fetched record, and we cannot set fields of a new sample either (ValueError).
		if sample not in fh_vcf_in.header.samples:
			new_rec = fh_vcf_out.header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=rec.alleles,
												   id=rec.id, qual=rec.qual, filter=rec.filter, info=rec.info)
			rec = new_rec

		target_sv = sv_class(rec, len_ratio_tol, ins_len_thr, del_len_thr, del_recip_overlap_thr, bnd_pos_tol)

		chrom = target_sv.chrom
		start = target_sv.start
		stop = target_sv.stop
		sv_id = target_sv.id
		svtype = target_sv.svtype
		svlen = target_sv.svlen
		chr2 = target_sv.chr2
		pos2 = target_sv.pos2

		if skip.skip_region(chrom, start, stop):
			count_skip_region += 1
			rec.info['SKIP_REGION'] = True
			fh_vcf_out.write(rec)
			continue

		if ((include_svtype != None and svtype not in include_svtype) or
			(exclude_svtype != None and svtype in exclude_svtype)):
			count_skip_svtype += 1
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

		read_supp_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()}
		read_supp_P_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### paternal
		read_supp_M_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### maternal
		read_supp_N_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### not phased

		if svlen != None and abs(svlen) > 10000:
			# we have two fetches here
			iter1 = fh_bam.fetch(chrom,
								 max(0, target_sv.start - buffer_length),
								 target_sv.start + buffer_length)
			iter2 = fh_bam.fetch(chrom,
								 max(0, target_sv.stop - buffer_length),
								 target_sv.stop + buffer_length)
		else:
			# we have one fetch here
			iter1 = fh_bam.fetch(chrom,
								 max(0, target_sv.start - buffer_length),
								 target_sv.stop + buffer_length)
			iter2 = iter(())
		for read in itertools.chain(iter1, iter2):
			if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
				locus_read, CG_supp, SA_supp = target_sv.sv_signature(read, mapping_quality_thr, buffer_length)
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
						raise NameError(f'problem with HP in read: {read.query_name}')
				else:
					read_supp_N_dict['locus_reads'].update([locus_read])
					read_supp_N_dict['CG_supp'].update([CG_supp])
					read_supp_N_dict['SA_supp'].update([SA_supp])

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
		GT, GQ, p_11, p_01, p_00, SQ = infer_gt_sv(DR_s, DV_s, p_err=p_err)
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
		print('count_skip_svtype:', count_skip_svtype)
	fh_bam.close()
	fh_vcf_in.close()
	fh_vcf_out.close()

