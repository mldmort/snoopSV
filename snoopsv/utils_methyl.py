from snoopsv.utils import get_ml

def methylation_signature(read, start, end, flanking_bp, verbose_debug):

	#print('I am here in the methylation signature function...')
	methyl_probs, num_bp = get_ml(read, start, end)

	return methyl_probs, num_bp
