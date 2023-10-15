import pysam
import subprocess
from snoopsv.utils import skip_region
from snoopsv.utils_sv import sv_class, sv_signature, infer_gt_sv, get_phased_gt
from snoopsv.utils_vcf import add_header_lines
from pathlib import Path

def GT_nonTR(vcf_in, vcf_out, contig, sample_bam_file, n_sec, i_sec, verbose=1):

	### genotyping setting
	mapping_quality_thr = 20
	region_buffer_length = 500
	SV_p_err = 0.01

	sample_bam_dict = {}
	if isinstance(sample_bam_file, str) or isinstance(sample_bam_file, Path):
		with open(sample_bam_file, 'r') as fh:
			for line in fh.readlines():
				line = line.strip().split()
				sample_bam_dict[line[0]] = line[1]
	elif isinstance(sample_bam_file, dict):
		sample_bam_dict = sample_bam_file
	else:
		raise ValueError(f'sample_bam_file variable does not have a good type: {type(sample_bam_file)}')

	sample2bamfh = {sample: pysam.AlignmentFile(bam, 'rb') for sample, bam in sample_bam_dict.items()}

	if contig:
		command = ('bcftools query -r '+contig+' -f %CHROM\\n '+vcf_in).split()
	else:
		command = ('bcftools query -f %CHROM\\n '+vcf_in).split()
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
	header_in = fh_vcf_in.header
	#for line in new_header_INFO+new_header_FORMAT:
	#	header_in.add_line(line)
	add_header_lines(header_in)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	count_skip_region = 0
	count_skip_sec = 0
	#count_skip_tr = 0
	for i_rec, rec in enumerate(fh_vcf_in.fetch(contig=contig)):
		if (i_rec < i_rec_start) or (i_rec >= i_rec_end):
			count_skip_sec += 1
			continue
		sv_id = rec.id
		svtype = rec.info['SVTYPE']
		svlen = rec.info['SVLEN']
		chrom = rec.chrom
		start = rec.start
		stop = rec.stop
		if svtype == 'TRA':
			chr2 = rec.info['CHR2']
		else:
			chr2 = chrom

		if skip_region(chrom, start, stop):
			count_skip_region += 1
			rec.info['SKIP_REGION'] = True
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
		target_sv = sv_class(rec)

		if svtype == 'INS':
			#pos_start = target_sv.start - 50
			#pos_stop = target_sv.stop + 49
			# don't need to pad INS calls when we don't intersect calls with TRs
			pos_start = target_sv.start
			pos_stop = target_sv.stop
		elif svtype != 'TRA':
			pos_start = target_sv.start
			pos_stop = target_sv.stop
		else:
			pos_start = target_sv.start
			pos_stop = target_sv.start+1

		for sample, bam_file in sample_bam_dict.items():
			#fh_bam = pysam.AlignmentFile(bam_file, 'rb')
			fh_bam = sample2bamfh[sample]
			read_supp_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()}
			read_supp_P_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### paternal
			read_supp_M_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### maternal
			read_supp_N_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### not phased
			read_supp_HF_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### HiFi reads
			read_supp_LF_dict = {'locus_reads':set(), 'CG_supp':set(), 'SA_supp':set()} ### LoFi reads

			for i_read, read in enumerate(fh_bam.fetch(chrom, max(0,pos_start-region_buffer_length), pos_stop+region_buffer_length)):
				if (not read.is_secondary) and (read.mapping_quality >= mapping_quality_thr):
					locus_read, CG_supp, SA_supp = sv_signature(read, target_sv)
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
							assert 0==1, 'problem with HP in read: ' + read.query_name
					else:
						read_supp_N_dict['locus_reads'].update([locus_read])
						read_supp_N_dict['CG_supp'].update([CG_supp])
						read_supp_N_dict['SA_supp'].update([SA_supp])
					if read.has_tag('HF'):
						HF = read.get_tag(tag='HF')
						if HF == 1:
							read_supp_HF_dict['locus_reads'].update([locus_read])
							read_supp_HF_dict['CG_supp'].update([CG_supp])
							read_supp_HF_dict['SA_supp'].update([SA_supp])
						elif HF == 0:
							read_supp_LF_dict['locus_reads'].update([locus_read])
							read_supp_LF_dict['CG_supp'].update([CG_supp])
							read_supp_LF_dict['SA_supp'].update([SA_supp])
						else:
							assert 0==1, 'problem with HF in read: ' + read.query_name
			#fh_bam.close()
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
			read_supp_HF_dict['locus_reads'] -= {''}
			read_supp_HF_dict['CG_supp'] -= {''}
			read_supp_HF_dict['SA_supp'] -= {''}
			read_supp_LF_dict['locus_reads'] -= {''}
			read_supp_LF_dict['CG_supp'] -= {''}
			read_supp_LF_dict['SA_supp'] -= {''}
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
			DV_s_HF = len(read_supp_HF_dict['CG_supp'] | read_supp_HF_dict['SA_supp'])
			DR_s_HF = len(read_supp_HF_dict['locus_reads']) - DV_s_HF
			assert DR_s_HF >= 0, 'problem with P DR/DV, DR: '+str(DR_s_HF)+', DV: '+str(DV_s_HF)+', sv_id: '+str(sv_id)
			DV_s_LF = len(read_supp_LF_dict['CG_supp'] | read_supp_LF_dict['SA_supp'])
			DR_s_LF = len(read_supp_LF_dict['locus_reads']) - DV_s_LF
			assert DR_s_LF >= 0, 'problem with P DR/DV, DR: '+str(DR_s_LF)+', DV: '+str(DV_s_LF)+', sv_id: '+str(sv_id)
			GT, GQ, p_11, p_01, p_00, SQ = infer_gt_sv(DR_s, DV_s, p_err=SV_p_err)
			GT_PH = get_phased_gt(GT, DV_s_P, DV_s_M)
			rec.samples[sample]['RV'] = DV_s
			rec.samples[sample]['RR'] = DR_s
			rec.samples[sample]['RV_P'] = DV_s_P
			rec.samples[sample]['RR_P'] = DR_s_P
			rec.samples[sample]['RV_M'] = DV_s_M
			rec.samples[sample]['RR_M'] = DR_s_M
			rec.samples[sample]['RV_N'] = DV_s_N
			rec.samples[sample]['RR_N'] = DR_s_N
			rec.samples[sample]['RV_HF'] = DV_s_HF
			rec.samples[sample]['RR_HF'] = DR_s_HF
			rec.samples[sample]['RV_LF'] = DV_s_LF
			rec.samples[sample]['RR_LF'] = DR_s_LF
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
		#print('count_skip_tr:', count_skip_tr)
	for sample, bam_fh in sample2bamfh.items():
		bam_fh.close()
	fh_vcf_in.close()
	fh_vcf_out.close()

