import pysam
import pandas as pd
import numpy as np
import pickle

def get_data(files, columns=None, header=None, index_col=None):

	data = pd.read_table(files[0], header=header, index_col=index_col, sep='\t', keep_default_na=False)
	for f in files[1:]:
		data = pd.concat([data, pd.read_table(f, header=None, index_col=None, sep='\t', keep_default_na=False)], ignore_index=True, axis=0)
	if columns != None:
		data.columns = columns
	return data

def get_data_lin(annot_in, cov_in):

	sample_cov_dict = {}
	with open(cov_in, 'r') as fh:
		for line in fh.readlines():
			sample = line.strip().split('\t')[0]
			cov_file = line.strip().split('\t')[1]
			sample_cov_dict[sample] = get_data([cov_file], header=0, index_col=0)

	main_chroms_list = ['chr'+str(x) for x in range(1,23)] + ['chrX', 'chrY']

	data_annot = pd.read_table(annot_in, sep='\t', header=0)
	#print(data_annot)
	data_annot = data_annot.loc[data_annot.chr2.isin(main_chroms_list)]
	#print(data_annot)

	if data_annot.empty:
		return pd.DataFrame()

	### data_lin will have columns_common columns plus these ones: rv, dp, dp_fr, mapq, sample
	columns_common = ['chrom', 'pos', 'end', 'id', 'svtype', 'chr2', 'end2']

	data_lin = pd.DataFrame()
	for index in data_annot.index.tolist():
		svtype = data_annot.loc[index, 'svtype']
		samples_list = str(data_annot.loc[index, 'samples']).split(',')
		mapq_list = str(data_annot.loc[index, 'mapq_all_samples']).split(',')
		dp_fr_list = str(data_annot.loc[index, 'dp_fr_samples']).split(',')
		dp_list = str(data_annot.loc[index, 'dp_samples']).split(',')
		rv_list = str(data_annot.loc[index, 'rv']).split(',')
		for i_sam, sample in enumerate(samples_list):
			mapq = mapq_list[i_sam]
			if svtype!='TRA':
				dp_fr = dp_fr_list[i_sam]
			else:
				dp_fr = '.'
			dp = dp_list[i_sam]
			rv = rv_list[i_sam]
			temp = pd.DataFrame(data_annot.loc[[index], columns_common])
			temp['rv'] = rv
			temp['dp'] = dp
			temp['dp_fr'] = dp_fr
			temp['mapq'] = mapq
			temp['sample'] = sample
			if svtype=='TRA':
				temp = pd.concat([temp, temp], axis=0, ignore_index=True)
				temp.loc[0, 'svtype'] = 'TRA_1'
				temp.loc[1, 'svtype'] = 'TRA_2'
				temp.loc[0, 'dp'] = temp.loc[0, 'dp'].split(':')[0]
				temp.loc[1, 'dp'] = temp.loc[1, 'dp'].split(':')[1]
			data_lin = pd.concat([data_lin, temp], axis=0, ignore_index=True)

	data_lin['rv'] = data_lin['rv'].astype(float)
	data_lin['dp'] = data_lin['dp'].astype(float)
	data_lin['mapq'] = data_lin['mapq'].astype(float)

	sub_data = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2')]
	if not sub_data.empty:
		### split dp_fr
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_1'] = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr'].str.split(pat=':', expand=True)[0]
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_2'] = data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr'].str.split(pat=':', expand=True)[1]
		### correct depth
		data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_cor'] = (data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_1'].astype(float) + data_lin.loc[(data_lin.svtype!='TRA_1') & (data_lin.svtype!='TRA_2'), 'dp_fr_2'].astype(float) + 1e-7) / 2.

	sub_data = data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2')]
	if not sub_data.empty:
		### correct depth
		data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2'), 'dp_cor'] = data_lin.loc[(data_lin.svtype=='TRA_1') | (data_lin.svtype=='TRA_2'), 'dp'] + 1e-7

	### compute DV/depth with corrected depth
	data_lin['rv/depth'] = data_lin.rv / data_lin.dp_cor

	### compute depth/cov
	for key in sample_cov_dict:
		for chrom in ['chr'+str(n) for n in range(1,23)]+['chrX', 'chrY']:
			data_lin.loc[(data_lin['sample']==key) & (data_lin.chrom==chrom) & (data_lin.svtype!='TRA_2'), 'cov'] = sample_cov_dict[key].loc[chrom, 'mean']
			data_lin.loc[(data_lin['sample']==key) & (data_lin.chr2==chrom) & (data_lin.svtype=='TRA_2'), 'cov'] = sample_cov_dict[key].loc[chrom, 'mean']
	data_lin['depth/cov'] = data_lin['dp_cor'] / data_lin['cov']

	return data_lin

def apply_ML_models(data_lin, model_labels, model_files, model_vars):

	for i_model, model_file in enumerate(model_files):

		model_var_list = model_vars[i_model]
		model_label = model_labels[i_model]
		print('working on model:', model_label)

		with open(model_file, 'rb') as fh:
			clf = pickle.load(fh)

		X_test = np.array(data_lin[model_var_list])
		predict_proba = clf.predict_proba(X_test)[:,1] # [:,0]: probablity of class 0; [:,1]: probablity of class 1

		data_lin[model_label] = predict_proba

	return data_lin

def write_vcf(vcf_in, data_lin, model_labels, vcf_out):

	new_headers = []
	for model_label in model_labels:
		new_headers.append('##FORMAT=<ID='+model_label+',Number=1,Type=String,Description="ML model score">')

	new_headers.append('##FORMAT=<ID=DP,Number=1,Type=String,Description="regional read depth">')
	new_headers.append('##FORMAT=<ID=COV,Number=1,Type=String,Description="average chromosome coverage">')
	new_headers.append('##FORMAT=<ID=MapQ,Number=1,Type=String,Description="average mapping quality of all intersecting reads">')
	new_headers.append('##FORMAT=<ID=DP_o_COV,Number=1,Type=String,Description="regional read depth over chrom coverage">')
	new_headers.append('##FORMAT=<ID=RV_o_DP,Number=1,Type=String,Description="read support over regional read depth">')
	new_headers.append('##FORMAT=<ID=DP_FR,Number=1,Type=String,Description="read depth of flanking regions">')

	fh_vcf_in = pysam.VariantFile(vcf_in)
	header_in = fh_vcf_in.header
	for line in new_headers:
		header_in.add_line(line)

	fh_vcf_out = pysam.VariantFile(vcf_out, mode='w', header=header_in)

	for i_rec, rec in enumerate(fh_vcf_in.fetch()):
		sv_id = rec.id
		svtype = rec.info['SVTYPE']

		if not data_lin.empty:
			if svtype != 'TRA':
				sub_data = data_lin.loc[data_lin.id==sv_id]
				for index in sub_data.index.tolist():
					sample = sub_data.loc[index, 'sample']
					dp = sub_data.loc[index, 'dp_cor']
					dp_fr = '|'.join(sub_data.loc[index, 'dp_fr'].split(':'))
					cov = sub_data.loc[index, 'cov']
					dp_o_cov = sub_data.loc[index, 'depth/cov']
					mapq = sub_data.loc[index, 'mapq']
					rv_o_dp = sub_data.loc[index, 'rv/depth']
					model_scores = {}
					for model in model_labels:
						model_scores[model] = sub_data.loc[index, model]

					rec.samples[sample]['DP'] = '{:.2f}'.format(dp)
					rec.samples[sample]['COV'] = '{:.2f}'.format(cov)
					rec.samples[sample]['MapQ'] = '{:.2f}'.format(mapq)
					rec.samples[sample]['DP_o_COV'] = '{:.2f}'.format(dp_o_cov)
					rec.samples[sample]['RV_o_DP'] = '{:.2f}'.format(rv_o_dp)
					rec.samples[sample]['DP_FR'] = str(dp_fr)
					for model, score in model_scores.items():
						rec.samples[sample][model] = '{:.2f}'.format(score)
			else:
				sub_data_1 = data_lin.loc[(data_lin.id==sv_id) & (data_lin.svtype=='TRA_1')]
				sub_data_2 = data_lin.loc[(data_lin.id==sv_id) & (data_lin.svtype=='TRA_2')]
				#print('sub_data_1:')
				#print(sub_data_1)
				#print('sub_data_2:')
				#print(sub_data_2)
				for index in sub_data_1.index.tolist():
					sample = sub_data_1.loc[index, 'sample']
					mapq = sub_data_1.loc[index, 'mapq']
					dp_1 = sub_data_1.loc[index, 'dp_cor']
					dp_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'dp_cor'].values[0]
					cov_1 = sub_data_1.loc[index, 'cov']
					cov_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'cov'].values[0]
					dp_o_cov_1 = sub_data_1.loc[index, 'depth/cov']
					dp_o_cov_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'depth/cov'].values[0]
					rv_o_dp_1 = sub_data_1.loc[index, 'rv/depth']
					rv_o_dp_2 = sub_data_2.loc[sub_data_2['sample']==sample, 'rv/depth'].values[0]
					model_scores_1 = {}
					model_scores_2 = {}
					for model in model_labels:
						model_scores_1[model] = sub_data_1.loc[index, model]
						model_scores_2[model] = sub_data_2.loc[sub_data_2['sample']==sample, model].values[0]

					rec.samples[sample]['DP'] = '{:.2f}'.format(dp_1)+'|'+'{:.2f}'.format(dp_2)
					rec.samples[sample]['COV'] = '{:.2f}'.format(cov_1)+'|'+'{:.2f}'.format(cov_2)
					rec.samples[sample]['MapQ'] = '{:.2f}'.format(mapq)
					rec.samples[sample]['DP_o_COV'] = '{:.2f}'.format(dp_o_cov_1)+'|'+'{:.2f}'.format(dp_o_cov_2)
					rec.samples[sample]['RV_o_DP'] = '{:.2f}'.format(rv_o_dp_1)+'|'+'{:.2f}'.format(rv_o_dp_2)
					for model, score_1 in model_scores_1.items():
						score_2 = model_scores_2[model]
						rec.samples[sample][model] = '{:.2f}'.format(score_1)+'|'+'{:.2f}'.format(score_2)
			
		fh_vcf_out.write(rec)

	fh_vcf_in.close()
	fh_vcf_out.close()
			
	
def SCORE_VCF(vcf_in, annot_in, cov_in, models, vcf_out):

	data_lin = get_data_lin(annot_in, cov_in)
	#print('data_lin:')
	#print(data_lin)

	model_labels = []
	model_files = []
	model_vars = []
	with open(models, 'r') as fh:
		for line in fh.readlines():
			line = line.strip().split('\t')
			model_labels.append(line[0])
			model_files.append(line[1])
			model_vars.append(line[2].split(','))
	#print('model_labels:', model_labels)
	#print('model_files:', model_files)
	#print('model_vars:', model_vars)
	if not data_lin.empty:
		data_lin = apply_ML_models(data_lin, model_labels, model_files, model_vars)
	#print('data_lin:')
	#print(data_lin)

	write_vcf(vcf_in, data_lin, model_labels, vcf_out)
