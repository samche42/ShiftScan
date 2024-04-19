#!/usr/bin/env python3

import pandas as pd
import numpy as np
from pathlib import Path
import argparse

parser = argparse.ArgumentParser()

#Add in input flags
parser.add_argument("-i", "--input_dir", help="Full path to input files")
parser.add_argument("-b", "--basename", help="Common name for files")

args = parser.parse_args()

#Get the experimental files from the path provided
files = Path(args.input_dir).glob('*.txt')  #generate list of all txt files
dfs = list() #Empty list of all files (converted to dfs)
for file in files:
    data = pd.read_csv(file, sep='\t',header=0)
    data['Origin of data'] = file.stem #Stem: Get filename without extension
    dfs.append(data)

raw_input_df = pd.concat(dfs, ignore_index=True)
raw_input_df.columns = raw_input_df.columns.str.replace("X.", "X (", regex = True)
raw_input_df.columns = [x+')' if 'X (' in x else x for x in raw_input_df.columns]

raw_input_df.to_csv(args.input_dir+'/'+args.basename+'_concatenated.txt', 
sep="\t", index=False)
