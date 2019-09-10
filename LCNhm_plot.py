"""
This script is quite simple; it takes data from the :data:`LCNhm_main.FolderName` ``LCNhm-results``'s subfolder,
and plot it. It shows on the left, the cell morphology and the recording sites, and on the right the membrane
potentials. 

To properly use this, you must type on the terminal ``python LCNhm_plot.py`` followed by the name of the experiment that
you want to be plotted (:data:`LCNhm_configurationfile.DIR_LOCATION). 

Eg.:

>>> python LCNhm_plot.py D20190730_T1030_First_Test
"""
import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from LCNhm_configurationfile import DIR_LOCATION, CELLPROP_MORPHOLOGY

# Input: Name of experiment to plot (folder inside ./LCNhm-results/)
Folder2plot = sys.argv[1]

# Names of the files that contains the Vmems, Recording positions, and Morphological data (SWC)
FileVmem = '%s/LCNhm-results/%s/Recordings_Vmem.txt'%(DIR_LOCATION,Folder2plot)
FilePos = '%s/LCNhm-results/%s/Recordings_Pos.txt'%(DIR_LOCATION,Folder2plot)
FileParam = '%s/LCNhm-results/%s/Parameters.txt'%(DIR_LOCATION,Folder2plot)
FileSWC = '%s/LCNhm-neurondata/%s.swc'%(DIR_LOCATION,CELLPROP_MORPHOLOGY)

# Get Vmems, Recording positions, and Morphological data (SWC)
Name = []
Vmem = []
Pos = []
xs, ys, zs, ds = [], [], [], []
# - Vmems
with open(FileVmem) as file:
	for ii, line in enumerate(file):
		if ii==0: 
			Time = map(float,line.split(' ')[1:])
		else: 
			Name.append( map(str,line.split(' ')[0]) )
			Vmem.append( map(float,line.split(' ')[1:]) )
# - Position
with open(FilePos) as file:
	for line in file:
		Pos.append( map(float,line.split(' ')[1:]) )
# - SWC
with open(FileSWC) as file:
	for line in file:
		x, y, z, d = map(float,line.split(' ')[2:6])
		xs.append(x)
		ys.append(y)
		zs.append(z)
		ds.append(d)

# Plot
f, (ax1, ax2, ax3) = plt.subplots(1, 3, gridspec_kw={'width_ratios': [1,2,1]}, figsize=(20,10))
ax1.scatter(xs, zs,color='k',s=ds)
ax1.axis('equal')
for ii in range(len(Pos)):
	color = np.random.rand(3)
	ax1.plot(Pos[ii][0]+20,Pos[ii][2],'<',color=color,markersize=20,alpha=0.6,markeredgecolor='k')
	ax2.plot(Time,Vmem[ii],color=color)
ax1.set_xlabel(u'x (\u03BCm)')
ax1.set_ylabel(u'y (\u03BCm)')
ax1.set_title('Recording places')
ax2.set_xlabel('Time (ms)')
ax2.set_ylabel('Vmem (mV)')
ax2.set_title('Membrane potential')
Parameters = ''
with open(FileParam) as file:
    for line in file:
        Parameters += line
ax3.text(0,0,Parameters)
ax3.axis('off')
plt.suptitle('Vmem from cell %s on simulation %s'%(CELLPROP_MORPHOLOGY,Folder2plot))
plt.savefig('%s/LCNhm-results/%s/Recordings_Plot.png'%(DIR_LOCATION,Folder2plot))
plt.show()