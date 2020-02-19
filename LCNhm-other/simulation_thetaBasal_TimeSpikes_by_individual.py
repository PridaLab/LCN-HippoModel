import os
import sys
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from LCNhm_configurationfile import DIR_LOCATION, SIMPROP_THETA_PERIOD, SIMPROP_SIM_TIME
DirLoc = np.copy(DIR_LOCATION)
tTheta = np.copy(SIMPROP_THETA_PERIOD)

Morphs = ['sup1','n409','n127','n128']#,'CA196f57-sup','cck-28-sup','PV-18-sup','CA-205R-F57-deep-mEC','CA211-LF-59-deep-AMG','CA211-RF-57-deep-mPFC','cck1-deep','PV20-deep','PV-63-deep-mEC','PV96-deep-mEC']
Intrs = [0, 1, 2, 4, 6, 7, 9, 12, 14, 15, 16, 19, 20, 21]
Syns = [0, 1, 3]

Spikes = {Morph: [ [ [] for Syn in Syns] for Intr in Intrs] for Morph in Morphs}
numCycs = []
with open(DIR_LOCATION+'/LCNhm-results/TimeSpikes_by_individual.txt') as file:
	for line in file:
		Intr = int(line.split(' ')[0]); iI = Intrs.index(Intr)
		Syn = int(line.split(' ')[1]); iS = Syns.index(Syn)
		Morph = line.split(' ')[2]
		if Morph in Morphs:
			if all(np.array(line[:-1].split(' ')[3:])==''): numCycs.append(len(line[:-1].split(' ')[3:]))
			spks = np.array(map(float,[ii for ii in line[:-1].split(' ')[3:] if ii!='']))
			Spikes[Morph][iI][iS].append(spks[spks>60.])
numCycs = np.mean(numCycs)
TimeTotal = numCycs * (SIMPROP_SIM_TIME-60) / 1000.

plt.figure(figsize=(15,40))
gs = gridspec.GridSpec(len(Intrs)*len(Syns),len(Morphs))
gs.update(wspace=0.025, hspace=0.05) # set the spacing between axes. 
k = 0
for Intr in Intrs:
	for Syn in Syns:
		for Morph in Morphs:
			iI = Intrs.index(Intr)
			iS = Syns.index(Syn)
			ax = plt.subplot(gs[k])
			ax.set_xticklabels([])
			ax.set_yticklabels([])
			Phases = np.mod( np.array(Spikes[Morph][iI][iS]) / tTheta * 360, 360)
			[yhist,xhist] = np.histogram(np.append(Phases,Phases+360), np.arange(0,720,15))
			plt.bar(xhist[:-1], yhist*0.99999/np.max(yhist), 15, color='k')
			plt.plot(np.arange(720),0.25*(1-np.cos(np.arange(720)/180.*np.pi)),color=[.5,.5,.5])
			plt.text(45,0.95,'%.2f Hz'%(len(Spikes[Morph][iI][iS][0])/TimeTotal),color='b')
			plt.ylim([0,1.2])
			if (iI==0) and (iS==0): ax.set_title(Morph)
			if (iI==len(Intrs)): ax.set_xticks(np.arange(0,721,180))
			if Morph==Morphs[0]: ax.set_ylabel('(%d,%d)'%(Intr,Syn))
			k+=1


plt.savefig(DIR_LOCATION+'/LCNhm-results/hist_TimeSpikes_by_individual_from60ms_2.png')
plt.show()