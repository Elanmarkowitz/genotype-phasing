import matplotlib
import json
from statistics import mean
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd

with open("hyperresults1.txt","r") as f:
	hyp1 = json.loads(f.read())

with open("hyperresults2.txt","r") as f:
	hyp2 = json.loads(f.read())

window = []
overlap = []
acc = []
for w in hyp1.keys():
	for o in hyp1[w].keys():
		m = mean([mean(hyp1[w][o].values()), mean(hyp2[w][o].values())])
		acc.append(m)
		window.append(int(w))
		overlap.append(int(o))

arr = pd.DataFrame([window, overlap, acc], ['window','overlap','acc'])
arr = arr.sort_values('acc', axis=1)

print(arr)

