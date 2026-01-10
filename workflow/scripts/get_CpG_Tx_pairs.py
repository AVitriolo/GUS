import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--input_dir_selected_CpGs', dest='input_dir_selected_CpGs', type=int, help='Add input_dir_selected_CpGs')
parser.add_argument('--output_path_pairs', dest='output_path_pairs', type=int, help='Add output_path_pairs')

args = parser.parse_args()

# load all files of selected features
selectedfol = args.input_dir_selected_CpGs
allfiles = os.listdir(selectedfol)

# take genes from last chunk of file name
genes = [x.split("_")[15] for x in allfiles]

cpgbygenetable = list() # empty list


# use idx of files list 
# loop over each file and each line
# paste gene name after selected feature

for idx in range(len(allfiles)):
	file = allfiles[idx]
	gene = genes[idx]
	with open(selectedfol+"/"+file,"r") as f:
		lines=f.readlines()
		for line in lines:
			line = line.strip()
			cpgbygenetable.append(line+"\t"+gene+"\n")

with open(args.output_path_pairs, "w") as f:
	f.writelines(cpgbygenetable)
