import numpy as np

def AlignmentScore(v, w, k_s_dict):
	if len(w) == 0:
		return [], []
	sigma = 1
	match = lambda x: 1 if x[0]==x[1] else -1
	Score_mat = np.zeros((len(v) + 1, len(w) + 1), dtype = int)
	### Backtrack convention:
	# 0: down
	# 1: right
	# 2: diag
	# 3: taxi from origin
	Backtrack = np.ones((len(v) + 1, len(w) + 1), dtype = int)
	Backtrack[0,:] = 1
	Backtrack[:,0] = 0
	for i in range(1, len(v) + 1):
		Score_mat[i, 0] = Score_mat[i - 1, 0] - sigma

	for i in range(1, len(v) + 1):
		for j in range(1, len(w) + 1):
			### options: list of options, each element: (score, direction)
			options = [(0, 3), (Score_mat[i - 1, j] - sigma, 0), (Score_mat[i, j - 1] - sigma, 1), (Score_mat[i - 1, j - 1] + match([v[i-1], w[j-1]]), 2)]
			Score_mat[i, j], Backtrack[i, j] = sorted(options, key=lambda x:x[0], reverse=True)[0]
	s_thr = k_s_dict[len(v)]
	#print('s_thr:', s_thr)
	score_last_row = Score_mat[-1,:]
	#print('score_last_row:', score_last_row)
	score_ind_list = []
	for ind, score in enumerate(score_last_row[1:]):
		if score >= s_thr:	
			score_ind_list.append((score, ind+1))
	#print('score_ind_list:', score_ind_list)

	if len(score_ind_list) == 0:
		return [], []

	cor_score_ind_list = []
	#print('cor_score_ind_list:')
	#print(cor_score_ind_list)
	ind_diff_thr = min(2, len(v)/2)
	score_ind_connected_list = [score_ind_list[0]]
	for i_score_ind_pre, score_ind_tuple in enumerate(score_ind_list[1:]):
		score, ind = score_ind_tuple
		score_pre, ind_pre = score_ind_list[i_score_ind_pre] # since we start from index 1 of score_ind_list
		if (ind - ind_pre) <= ind_diff_thr:
			score_ind_connected_list.append((score, ind))
		else:
			cor_score_ind_list.append(sorted(score_ind_connected_list, key=lambda x:x[0], reverse=True)[0])
			score_ind_connected_list = [(score, ind)]
	if len(score_ind_connected_list) > 0:
		cor_score_ind_list.append(sorted(score_ind_connected_list, key=lambda x:x[0], reverse=True)[0])
	#print('cor_score_ind_list:')
	#print(cor_score_ind_list)
	return cor_score_ind_list, score_ind_list

