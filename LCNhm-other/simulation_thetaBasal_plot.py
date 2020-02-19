import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from LCNhm_configurationfile import DIR_LOCATION
DirLoc = np.copy(DIR_LOCATION)


def compute_prefer_phase(spktimes):
	
	tTheta = 166.

	spktimes = np.array(spktimes)
	spktimes = np.concatenate(spktimes).ravel()
	spktimes = spktimes[~np.isnan(spktimes)]
	if len(spktimes)>0:
		yhist, xhist = np.histogram((spktimes%tTheta)*360./tTheta,np.arange(0,360,5))
		imax = np.argmax(yhist)
		return xhist[imax]
	else:
		return np.nan






MorphNames = ['sup1','n409','n127','n128','cck-28-sup','PV-18-sup','CA170-F59-deep-AMG','I040913C2','PV20-deep','CA196f57-sup','I230114B14C4','PV-63-deep-mEC','CA-205R-F57-deep-mEC','PV75-deep-mEC','CA211-LF-59-deep-AMG','PV96-deep-mEC','CA211-RF-57-deep-mPFC','cck1-deep','PV-12-deep-mPFC']

SearchText = 'thetaBasal_test'

DirList = [Dir for Dir in os.listdir('./LCNhm-results/%s/'%SearchText) if SearchText in Dir]
Syns = []
Morphs = []
Intrs = []
Spks = [ [ [ [] for ii in range(4) ] for jj in range(19) ] for kk in range(21)]
Phases = [ [ [ [] for ii in range(4) ] for jj in range(19) ] for kk in range(21)]
Vs = [ [ [ [] for ii in range(4) ] for jj in range(19) ] for kk in range(21)]
VList = []

for Folder2plot in DirList:

	# Names of the files that contains the Vmems, Recording positions, and Morphological data (SWC)
	#FileVmem = '%s/LCNhm-results/%s/%s/Recordings_Vmem.txt'%(DirLoc,SearchText,Folder2plot)
	#FileParam = '%s/LCNhm-results/%s/%s/Parameters.txt'%(DirLoc,SearchText,Folder2plot)
	FileParam = '%s/LCNhm-results/%s/%s/Parameters.txt'%(DirLoc,SearchText,Folder2plot)
	FileSpks = '%s/LCNhm-results/%s/%s/TimeSpikes.txt'%(DirLoc,SearchText,Folder2plot)
	FileVmem = '%s/LCNhm-results/%s/%s/Recordings_Vmem.txt'%(DirLoc,SearchText,Folder2plot)

	if os.path.isfile(FileSpks):
		with open(FileParam) as file:
			for line in file:
				VarName = line.split(' ')[0]
				VarVals = line[:-1].split(' ')[1:]
				if len(VarVals)==1:
					exec("%s = %s" % (VarName,VarVals))
				else:
					try:
						exec("%s = %s" % (VarName,"["+", ".join(VarVals)+"]"))
					except:
						''

		if 'CELLPROP_SYNAPTIC' in globals():
			Intr = int(CELLPROP_INDIVIDUAL[0])
			Morph = [iMorph for iMorph in range(len(MorphNames)) if MorphNames[iMorph]==CELLPROP_MORPHOLOGY[0]][0]
			Syn = int(CELLPROP_SYNAPTIC[0])
			DT = float(SIMPROP_DT[0])

			# Get spikes
			'''
			Time = []
			V = []
			if os.path.isfile(FileVmem):
				with open(FileVmem) as file:
					for ii, line in enumerate(file):
						if ii==0: 
							Time = map(float,line.split(' ')[1:])
						else: 
							V.append( map(float,line.split(' ')[1:]) )
			'''

			# Get spikes
			Spk = np.array([])
			with open(FileSpks) as file:
				for line in file:
					Spk = np.append(Spk, map(float,line.split(' ')))

			Syns.append(Syn)
			Morphs.append(Morph)
			Intrs.append(Intr)
			#Vs[Intr][Morph][Syn].append(V)
			#VList.append(V)
			Spks[Intr][Morph][Syn].append(Spk)

# Get preferred phase
for iI in range(np.shape(Spks)[0]):
	for iM in range(np.shape(Spks)[1]):
		for iS in range(np.shape(Spks)[2]):
			Phases[iI][iM][iS] = compute_prefer_phase(Spks[iI][iM][iS]) if (len(Spks[iI][iM][iS])>0) else np.nan

Syns = np.array(Syns)
Morphs = np.array(Morphs)
Intrs = np.array(Intrs)
Phases = np.array(Phases)
Spks = np.array(Spks)
#Vs = np.array(Vs)
#VList = np.array(VList)
#Time = np.array(Time)









# Plot preferred phases
fig, axs = plt.subplots(1,4,figsize=(20,5))
phases2plot = (Phases + 65)%360 
for ii in range(4):
	img = axs[ii].imshow(np.transpose(phases2plot[:,:,ii]), vmin=0, vmax=360, cmap='coolwarm')
	#fig.colorbar(img)
	axs[ii].set_xticks(range(22))
	axs[ii].set_xticklabels(np.unique(Intrs))
	if ii==0: 
		axs[ii].set_yticks(range(len(MorphNames)))
		axs[ii].set_yticklabels(MorphNames)
	else:
		axs[ii].set_yticks([])
	axs[ii].set_title('Set synpases #%d'%ii)
plt.savefig('LCNhm-results/%s_phases.png'%SearchText)
plt.show()


# Plot Vmems
for morph2plot in range(7)+range(8,19):
	intrs2plot = [0,1,2,5,9,10,12,16,20,21]
	plt.figure(figsize=(25,10))
	plt.plot(Time,np.transpose(np.squeeze(VList[(Morphs==morph2plot) * np.isin(Intrs,intrs2plot)])))
	plt.plot(Time,-75+10*(1-np.cos(2*np.pi*Time/166.)),'k')
	plt.title(MorphNames[morph2plot])
	plt.xlabel('Time (ms)')
	plt.ylabel('Vmem (mV)')
	plt.savefig('LCNhm-results/%s_Vmem_Morph_%s.png'%(SearchText,MorphNames[morph2plot]))
plt.close()
plt.show()







