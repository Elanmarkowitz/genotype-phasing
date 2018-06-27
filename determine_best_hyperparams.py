import phaser as ph
import subprocess
import sys
import json

tester_script = 'calculate_switch_accuracy.R'

accuracies = {}
for window in [31, 35,39]:
	accuracies[window] = {}
	for overlap in [3, 5, 7]:
		if overlap*2 >= window:
			continue
		accuracies[window][overlap] = {}
		for num_seeds in [5,10]:
			ph.run_predictions(
				"example_data_2.txt", 
				"res2_{}_{}_{}".format(window, overlap, num_seeds), 
				window, 
				overlap, 
				num_seeds)
			run_check_cmd = ['Rscript', tester_script, "res2_{}_{}_{}".format(window, overlap, num_seeds), "example_data_2_sol.txt"]
			try:
				x = subprocess.check_output(run_check_cmd, universal_newlines=True)
			except:
				x=-1
				print("subprocess failed")
			try:
				accuracies[window][overlap][num_seeds] = float(x)
			except:
				print("error")
			f = open('hyperresults2.txt', 'w')
			f.write(json.dumps(accuracies))
			f.close()

accuracies = {}
for window in [31, 35,39]:
	accuracies[window] = {}
	for overlap in [3, 5, 7]:
		if overlap*2 >= window:
			continue
		accuracies[window][overlap] = {}
		for num_seeds in [5,10]:
			ph.run_predictions(
				"example_data_1.txt", 
				"res1_{}_{}_{}".format(window, overlap, num_seeds), 
				window, 
				overlap, 
				num_seeds)
			run_check_cmd = ['Rscript', tester_script, "res1_{}_{}_{}".format(window, overlap, num_seeds), "example_data_1_sol.txt"]
			x = subprocess.check_output(run_check_cmd, universal_newlines=True)
			accuracies[window][overlap][num_seeds] = float(x)
			f = open('hyperresults1.txt', 'w')
			f.write(json.dumps(accuracies))
			f.close()


