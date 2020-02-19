import os
import sys
import numpy as np
import matplotlib.pyplot as plt

DIR_LOCATION = '/home/andrea/Projects/HippoModel/LCNhippomodel'
FileSWCs = ['CA-205R-F57-deep-mEC','CA170-F59-deep-AMG','CA211-LF-59-deep-AMG','CA211-RF-57-deep-mPFC','cck1-deep','PV20-deep','PV-12-deep-mPFC','PV-63-deep-mEC','PV75-deep-mEC','PV96-deep-mEC','I230114B14C4','CA196f57-sup','cck-28-sup','PV-18-sup']

for morph in FileSWCs:
	FileSWC = '%s/LCNhm-neurondata/%s.swc'%(DIR_LOCATION,morph)

	with open(FileSWC, 'r') as input_file, open('tmp.txt', 'w') as output_file:
	    for line in input_file:
	        ID, t, x, y, z, d, f = map(float,line.split(' '))
	        output_file.write(' '.join(map(str,[int(ID), int(t), x, z, y, d, int(f)]))+'\n')
	os.rename('tmp.txt',FileSWC)