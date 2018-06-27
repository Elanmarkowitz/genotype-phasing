import os
import sys
import numpy as np
import collections
import pandas as pd
import random
import itertools

from phasertools import *
#input_file = sys.argv[1]
#output_file = sys.argv[2]

def get_haplotypes_from_seed(g, seed):
	# assign first het to (0,1) and then according to seed, and then random
	h1, h2 = [], []
	first = True
	i = 0
	for s in g:
		if s == 0:
			h1.append(0)
			h2.append(0)
		elif s == 2:
			h1.append(1)
			h2.append(1)
		elif s == 1:
			if first:
				h1.append(0)
				h2.append(1)
				first = False
				continue
			if i < len(seed) and seed[i] == 0:
				h1.append(0)
				h2.append(1)
			elif i < len(seed) and seed[i] == 1:
				h1.append(1)
				h2.append(0)
			else: 
				print("random")
				if random.getrandbits(1):
					h1.append(0)
					h2.append(1)
				else:
					h1.append(1)
					h2.append(0)
				i += 1
		else:
			return -1
	return h1, h2

def get_seeds(genotypes, num_seeds):
	least_ambiguous_g = []
	count_1_min = 10000000000
	for g in genotypes:
		count_1 = collections.Counter(g)[1]
		if count_1 < count_1_min:
			least_ambiguous_g = g
			count_1_min = count_1
	seed_Hs = []
	seed_het_indexes = []
	if count_1_min < 2:
		seed_het_activation = [(0)]
	elif count_1_min < num_seeds:
		seed_het_activation = itertools.product([0, 1], repeat=count_1_min - 1)
	else:
		seed_het_activation = itertools.product([0, 1], repeat=num_seeds)
	for seed in seed_het_activation:
			h1, h2 = get_haplotypes_from_seed(least_ambiguous_g, seed)
			seed_Hs.append((h1,h2))
	return seed_Hs

def run_clark(genotypes, H):
	predictions = [[]]*len(genotypes)*2
	predicted = set()
	for i, g in enumerate(genotypes):
		if is_unambiguous(g):
			h1, h2 = get_unambiguous_haplotype(g)
			H.add(tuple(h1))
			H.add(tuple(h2))
			predictions[2*i] = h1
			predictions[2*i + 1] = h2
			predicted.add(i)
	# add all cases that are fully explained by H
	for i, g in enumerate(genotypes):
		if i in predicted:
			continue
		for h in H:
			if is_compatible_haplotype(h, g):
				h2 = corresponding_haplotype(h, g)
				if tuple(h2) in H:
					predictions[2*i] = h
					predictions[2*i + 1] = h2
					predicted.add(i)
					break
	# add all cases that are half explained by H
	# and add corresponding h to H
	# continue until no more can be done
	while(len(predicted) < len(genotypes)):
		new_h = True
		while new_h:
			new_h = False
			for i, g in enumerate(genotypes):
				if i in predicted:
					continue
				for h in H:
					if is_compatible_haplotype(h, g):
						h2 = corresponding_haplotype(h, g)
						H.add(tuple(h2))
						predictions[2*i] = h
						predictions[2*i + 1] = h2
						predicted.add(i)
						new_h = True
						break
			if not new_h:
				for i, g in enumerate(genotypes):
					if i not in predicted:
						h1, h2 = get_rand_compatible_haplotypes(g)
						predictions[2*i] = h1
						predictions[2*i + 1] = h2
						predicted.add(i)
						H.add(tuple(h1))
						H.add(tuple(h2))
						new_h = True
						break
	return predictions, H

def predict_subset(genotypes, num_seeds):
	# implementation of clark's algorithm, run multiple times if ambiguous start
	predictions = [[]]*len(genotypes)*2
	H = set()
	min_H = len(genotypes)*2
	min_preds = predictions
	# add all unambiguous cases
	for seed_haplotypes in get_seeds(genotypes, num_seeds):
		H = set()
		H.add(tuple(seed_haplotypes[0]))
		H.add(tuple(seed_haplotypes[1]))
		predictions, H = run_clark(genotypes, H)
		if len(H) < min_H:
			min_H = len(H)
			min_preds = predictions
	#print("haplotypes", min_H)
	return np.array(min_preds)

# def merge(haps_1, haps_2, window, overlap):
#     n_rows = haps_1.shape[0]
#     n_columns = haps_1.shape[1] + haps_2.shape[1] - overlap
#     to_merge = haps_2[:, overlap:window]
#     if len(to_merge.shape) == 1:
#         to_merge = to_merge.reshape(n_rows, 1)
#     merged = np.zeros((n_rows, n_columns), dtype = int)
#     for i in range(int(n_rows/2)):
#         if haps_1[2 * i, haps_1.shape[1] - 2] == haps_2[2 * i, 0] and haps_1[2 * i, haps_1.shape[1] - 1] == haps_2[2 * i, 1]:
#             merged[2 * i] = np.hstack((haps_1[2 * i], to_merge[2 * i]))
#             merged[2 * i + 1] = np.hstack((haps_1[2 * i + 1], to_merge[2 * i + 1]))
#         else:
#             merged[2 * i] = np.hstack((haps_1[2 * i], to_merge[2 * i + 1]))
#             merged[2 * i + 1] = np.hstack((haps_1[2 * i + 1], to_merge[2 * i]))
#     return merged

def merge(hap1, hap2, window, overlap):
	new_hap = np.empty((len(hap1), len(hap1[0])+len(hap2[0])-overlap))
	overlap_start_pos = len(hap1[0]) - overlap
	#print(hap1.shape, hap2.shape)
	for i in range(len(hap1)//2):
		for rel_pos	in range(overlap + 1):
			abs_pos = rel_pos + overlap_start_pos
			if (rel_pos == overlap):
				# what to do for now het site in overlap
				new_hap[2*i] = np.concatenate((hap1[2*i], hap2[2*i,overlap:]))
				new_hap[2*i + 1] = np.concatenate((hap1[2*i + 1], hap2[2*i + 1,overlap:]))
				#print('unknown alignment')
				break
			if (hap1[2*i][abs_pos] != hap1[2*i+1][abs_pos]): #het site
				if hap1[2*i][abs_pos] == hap2[2*i][rel_pos]:
					#print('aligned')
					new_hap[2*i] = np.concatenate((hap1[2*i], hap2[2*i,overlap:]))
					new_hap[2*i + 1] = np.concatenate((hap1[2*i + 1], hap2[2*i + 1,overlap:]))
				else:
					#print('unaligned')
					new_hap[2*i] = np.concatenate((hap1[2*i], hap2[2*i+1,overlap:]))
					new_hap[2*i + 1] = np.concatenate((hap1[2*i + 1], hap2[2*i,overlap:]))
				break
	return np.int_(new_hap)

def predict_haplotypes(data, window, overlap, num_seeds):
	subsets = subset_data(data, window, overlap)
	for i, d in enumerate(subsets):
		if i == 0:
			preds = predict_subset(d, num_seeds)
			#print(preds.shape)
		else:
			pred = predict_subset(d, num_seeds)
			#print(pred.shape)
			preds = merge(preds, pred, window, overlap)
		# Progress bar
		sys.stdout.write("\rProgress: [" + "#"*(30*(i+1)//len(subsets)) + " "*(30 - 30*(i+1)//len(subsets)) + "]")
	return preds


def run_predictions(file_in, pred_out, window, overlap, num_seeds):
	data = load_data(file_in)
	res = predict_haplotypes(data, window, overlap, num_seeds)
	pd.DataFrame(res.transpose()).to_csv(pred_out, index=False, header=False, sep = " ")

def main():
	print("Running algorithm.")
	file_in = sys.argv[1]
	file_out = sys.argv[2]
	optimal_window = 35
	optimal_overlap = 5
	optimal_num_seeds = 10
	run_predictions(file_in, file_out, optimal_window, optimal_overlap, optimal_num_seeds)

if __name__ == "__main__":
    main()

# pd.DataFrame(res.transpose()).to_csv('testres.txt', index=False, header=False, sep = " ")