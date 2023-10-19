from pathlib import Path
import sys
from snoopsv.snoop_sv import GT_nonTR
import tempfile
import pysam

DATA_DIR = Path(__file__).parent.joinpath('../data')
BAM_INS_DEL_INV_DUP = str(DATA_DIR.joinpath('reads_INS_DEL_INV_DUP.bam'))
BAM_DUP_LARGE = str(DATA_DIR.joinpath('reads_DUP_large.bam'))
VCF_DEL = str(DATA_DIR.joinpath('Ecoli_simulated_DEL.sorted.vcf.gz'))
VCF_INS = str(DATA_DIR.joinpath('Ecoli_simulated_INS.sorted.vcf.gz'))
VCF_DUP = str(DATA_DIR.joinpath('Ecoli_simulated_DUP.sorted.vcf.gz'))
VCF_INV = str(DATA_DIR.joinpath('Ecoli_simulated_INV.sorted.vcf.gz'))
VCF_DUP_LARGE= str(DATA_DIR.joinpath('Ecoli_simulated_DUP_large.sorted.vcf.gz'))

skip_bed = None
mapping_quality_thr = 20
buffer_length = 500
p_err = 0.01
len_ratio_tol = 0.25
ins_len_thr = 20
del_len_thr = 20
del_recip_overlap_thr = 0.3
extra_params = (skip_bed, mapping_quality_thr, buffer_length, p_err, len_ratio_tol, ins_len_thr, del_len_thr, del_recip_overlap_thr)

def test_sv_del():
	temp_dir = tempfile.TemporaryDirectory()
	vcf_in = VCF_DEL
	vcf_out = temp_dir.name + '/' + 'test_vcf_out.vcf'
	contig = None
	sample = 'SAMPLE'
	bam = BAM_INS_DEL_INV_DUP
	n_sec = 1
	i_sec = 0
	GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, *extra_params, verbose=1)

	vcf_fh = pysam.VariantFile(vcf_out)
	for record in vcf_fh.fetch():
		if record.id == 'Sniffles2.DEL.50S0':
			assert record.samples["SAMPLE"]["RR"] == 0
			assert record.samples["SAMPLE"]["RV"] == 9
		elif record.id == 'Sniffles2.DEL.54S0':
			assert record.samples["SAMPLE"]["RR"] == 3 #TODO: this should be looked at. are all reads that are not
			# supporting the variant supporting the reference actually? or they can be non-informative?
			assert record.samples["SAMPLE"]["RV"] == 3
		elif record.id == 'Sniffles2.DEL.5AS0':
			assert record.samples["SAMPLE"]["RR"] == 2
			assert record.samples["SAMPLE"]["RV"] == 5
		elif record.id == 'Sniffles2.DEL.61S0':
			assert record.samples["SAMPLE"]["RR"] == 0
			assert record.samples["SAMPLE"]["RV"] == 6
		else:
			assert 1 == 0
		#print(f'id: {record.id}, RR: {record.samples["SAMPLE"]["RR"]}, RV: {record.samples["SAMPLE"]["RV"]}' )
	vcf_fh.close()

def test_sv_ins():
	temp_dir = tempfile.TemporaryDirectory()
	vcf_in = VCF_INS
	vcf_out = temp_dir.name + '/' + 'test_vcf_out.vcf'
	contig = None
	sample = 'SAMPLE'
	bam = BAM_INS_DEL_INV_DUP
	n_sec = 1
	i_sec = 0
	GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, *extra_params, verbose=1)

	vcf_fh = pysam.VariantFile(vcf_out)
	for record in vcf_fh.fetch():
		if record.id == 'Sniffles2.INS.BS0':
			assert record.samples["SAMPLE"]["RR"] == 0
			assert record.samples["SAMPLE"]["RV"] == 8
		elif record.id == 'Sniffles2.INS.2AS0':
			assert record.samples["SAMPLE"]["RR"] == 4 #TODO: this should be looked at. are all reads that are not
			# supporting the variant supporting the reference actually? or they can be non-informative?
			assert record.samples["SAMPLE"]["RV"] == 6 #TODO: this should be looked at. check CIGAR and supp
			# evidance separately.
		elif record.id == 'Sniffles2.INS.2CS0':
			assert record.samples["SAMPLE"]["RR"] == 1 #TODO: this should be looked at. are all reads that are not
			# supporting the variant supporting the reference actually? or they can be non-informative?
			assert record.samples["SAMPLE"]["RV"] == 4
		elif record.id == 'Sniffles2.INS.33S0':
			assert record.samples["SAMPLE"]["RR"] == 3 #TODO: this should be looked at. are all reads that are not
			# supporting the variant supporting the reference actually? or they can be non-informative?
			assert record.samples["SAMPLE"]["RV"] == 2
		else:
			assert 1 == 0
		#print(f'id: {record.id}, RR: {record.samples["SAMPLE"]["RR"]}, RV: {record.samples["SAMPLE"]["RV"]}' )
	vcf_fh.close()

def test_sv_dup():
	temp_dir = tempfile.TemporaryDirectory()
	vcf_in = VCF_INV
	vcf_out = temp_dir.name + '/' + 'test_vcf_out.vcf'
	contig = None
	sample = 'SAMPLE'
	bam = BAM_INS_DEL_INV_DUP
	n_sec = 1
	i_sec = 0
	GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, *extra_params, verbose=1)

	vcf_fh = pysam.VariantFile(vcf_out)
	print(vcf_fh)
	for record in vcf_fh.fetch():
		if record.id == 'Sniffles2.INV.81S0':
			assert record.samples["SAMPLE"]["RR"] == 1 #TODO: is there more RR according to sniffles?
			assert record.samples["SAMPLE"]["RV"] == 3
		elif record.id == 'Sniffles2.INV.83S0':
			assert record.samples["SAMPLE"]["RR"] == 3
			assert record.samples["SAMPLE"]["RV"] == 2 # can we squeese one more RV from a non-complete signal?
		elif record.id == 'Sniffles2.INV.86S0':
			assert record.samples["SAMPLE"]["RR"] == 2 #TODO: can we squeese one more RV from a non-complete signal?
			assert record.samples["SAMPLE"]["RV"] == 8
		elif record.id == 'Sniffles2.INV.92S0':
			assert record.samples["SAMPLE"]["RR"] == 1
			assert record.samples["SAMPLE"]["RV"] == 5 #TODO: where does one of them come from?
		else:
			assert 1 == 0
		#print(f'id: {record.id}, RR: {record.samples["SAMPLE"]["RR"]}, RV: {record.samples["SAMPLE"]["RV"]}' )
	vcf_fh.close()

def test_sv_inv():
	temp_dir = tempfile.TemporaryDirectory()
	vcf_in = VCF_INV
	vcf_out = temp_dir.name + '/' + 'test_vcf_out.vcf'
	contig = None
	sample = 'SAMPLE'
	bam = BAM_INS_DEL_INV_DUP
	n_sec = 1
	i_sec = 0
	GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, *extra_params, verbose=1)

	vcf_fh = pysam.VariantFile(vcf_out)
	for record in vcf_fh.fetch():
		if record.id == 'Sniffles2.INV.81S0':
			assert record.samples["SAMPLE"]["RR"] == 1
			assert record.samples["SAMPLE"]["RV"] == 3
		elif record.id == 'Sniffles2.INV.83S0':
			assert record.samples["SAMPLE"]["RR"] == 3
			assert record.samples["SAMPLE"]["RV"] == 2
		elif record.id == 'Sniffles2.INV.86S0':
			assert record.samples["SAMPLE"]["RR"] == 2
			assert record.samples["SAMPLE"]["RV"] == 8
		elif record.id == 'Sniffles2.INV.92S0':
			assert record.samples["SAMPLE"]["RR"] == 1
			assert record.samples["SAMPLE"]["RV"] == 5
		else:
			assert 1 == 0
		#print(f'id: {record.id}, RR: {record.samples["SAMPLE"]["RR"]}, RV: {record.samples["SAMPLE"]["RV"]}' )
	vcf_fh.close()

def test_sv_dup_large():
	temp_dir = tempfile.TemporaryDirectory()
	vcf_in = VCF_DUP_LARGE
	vcf_out = temp_dir.name + '/' + 'test_vcf_out.vcf'
	contig = None
	sample = 'SAMPLE'
	bam = BAM_DUP_LARGE
	n_sec = 1
	i_sec = 0
	GT_nonTR(vcf_in, vcf_out, contig, sample, bam, n_sec, i_sec, *extra_params, verbose=1)

	vcf_fh = pysam.VariantFile(vcf_out)
	for record in vcf_fh.fetch():
		if record.id == 'Sniffles2.DUP.2S0':
			assert record.samples["SAMPLE"]["RR"] == 20
			assert record.samples["SAMPLE"]["RV"] == 6
		elif record.id == 'Sniffles2.DUP.3S0':
			assert record.samples["SAMPLE"]["RR"] == 24
			assert record.samples["SAMPLE"]["RV"] == 5
		elif record.id == 'Sniffles2.DUP.4S0':
			assert record.samples["SAMPLE"]["RR"] == 23
			assert record.samples["SAMPLE"]["RV"] == 12
		else:
			assert 1 == 0
		#print(f'id: {record.id}, RR: {record.samples["SAMPLE"]["RR"]}, RV: {record.samples["SAMPLE"]["RV"]}' )
	vcf_fh.close()

if __name__ == '__main__':
	#sys.exit(test_sv_del())
	#sys.exit(test_sv_ins())
	#sys.exit(test_sv_dup())
	sys.exit(test_sv_inv())
	#sys.exit(test_sv_dup_large())

