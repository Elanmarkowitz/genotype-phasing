import os
import sys
import numpy as np
import collections
import pandas as pd
import random
import itertools

def load_data(file_in):
	f = open(file_in, "r")
	line = f.readline()[:-1]
	data = []
	while(line):
		line = line.split(" ")
		line = [int(x) for x in line]
		data.append(line)
		line = f.readline()[:-1]
	return np.transpose(data)

def is_unambiguous(g):
	count = collections.Counter(g)[1]
	return count < 2

def contains_unambiguous(genotypes):
	for g in genotypes:
		if is_unambiguous(g):
			return True
	return False

def subset_data_by_unambiguous(data):
	if contains_unambiguous(data):
		return [data]
	numsites = data.shape[1]
	split_data = np.split(data, [numsites//2], 1)
	return subset_data_by_unambiguous(split_data[0]) + subset_data_by_unambiguous(split_data[1])
		# could be possibly improved by greedily looking for the largest subsets that contain unambiguous samples

def subset_data(data, window, overlap):
	step_size = window - overlap
	subsets = []
	cur_point = 0
	while cur_point	+ overlap < data.shape[1]:
		cur_set = np.split(data, [cur_point, cur_point + window], axis = 1)
		subsets.append(cur_set[1])
		cur_point += step_size
	return subsets

def get_unambiguous_haplotype(g):
	h1, h2 = [], []
	if not is_unambiguous(g):
		return -1
	for s in g:
		if s == 0:
			h1.append(0)
			h2.append(0)
		elif s == 2:
			h1.append(1)
			h2.append(1)
		elif s == 1:
			h1.append(0)
			h2.append(1)

		else:
			return -1
	return h1, h2

def corresponding_haplotype(h, g):
	h2 = []
	for i, s in enumerate(g):
		if s == 0:
			h2.append(0)
		elif s == 1:
			h2.append(0 if h[i] == 1 else 1)
		elif s == 2:
			h2.append(1)
		else:
			return -1
	return h2

def is_compatible_haplotype(hap, gen):
	for h, g in zip(hap, gen):
		if (h == 1 and g == 0) or (h == 0 and g == 2):
			return False
	return True

def get_rand_compatible_haplotypes(g):
	h1, h2 = [], []
	for s in g:
		if s == 0:
			h1.append(0)
			h2.append(0)
		elif s == 2:
			h1.append(1)
			h2.append(1)
		elif s == 1:
			if random.getrandbits(1):
				h1.append(0)
				h2.append(1)
			else:
				h1.append(1)
				h2.append(0)
		else:
			return -1
	return h1, h2