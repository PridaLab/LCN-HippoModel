import os
import sys
import numpy as np
from scipy import stats
import statsmodels.api as sm
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from sklearn.linear_model import LinearRegression
from LCNhm_configurationfile import DIR_LOCATION, SIMPROP_THETA_PERIOD, SIMPROP_SIM_TIME
DirLoc = np.copy(DIR_LOCATION)
tTheta = np.copy(SIMPROP_THETA_PERIOD)


def compute_prefer_phase(spktimes):
	
	tTheta = 166.

	spktimes = np.array(spktimes)
	spktimes = spktimes[~np.isnan(spktimes)]
	yhist, xhist = np.histogram((spktimes%tTheta)*360./tTheta,np.arange(0,360,20))
	imax = np.argmax(yhist)

	return xhist[imax]

def compute_pval(x,y):

	'''
	_, pval = stats.ttest_rel( group1, group2 )
	_, pval = stats.ttest_ind( x, y )
	'''
	x = np.array(x).reshape(-1,1);
	model = LinearRegression().fit(x,y)
	x2 = sm.add_constant(x)
	stat = sm.OLS(y, x2)
	stat = stat.fit()
	pval = stat.pvalues[1]
	return pval, stat.params[0], stat.params[1], stat.rsquared


# 		==================================
# 			DEFINING FACTORS 			 
# 		==================================

# Intrinsic factor values for selected Intrs in (0, 1, 2, 4, 6, 7, 9, 12, 14, 15, 16, 19, 20, 21)
# 				fNa,   fA,    fAHPs, fC,    fKDR,  fM,    fCa,   fHCN,  fL,    Ra
IntrFactors = [[6.458, 0.495, 3.811, 1.796, 2.500, 1.232, 1.800, 1.000, 1.000, 20.000 ], 
              [ 6.820, 0.500, 12.389, 2.500, 2.500, 25.000, 1.800, 1.000, 1.000, 20.000 ], 
              [ 7.500, 0.500, 27.138, 2.500, 2.500, 10.000, 1.800, 1.000, 1.000, 20.000 ], 
              [ 15.574, 0.500, 0.903, 26.605, 9.315, 25.000, 1.800, 25.000, 0.500, 20.000 ], 
              [ 18.500, 1.357, 0.162, 13.038, 25.383, 1.926, 4.211, 25.843, 2.490, 8.820 ], 
              [ 18.500, 0.700, 12.389, 7.919, 2.500, 25.000, 2.500, 1.000, 1.000, 20.000 ], 
              [ 21.088, 1.357, 0.505, 26.029, 26.862, 1.926, 29.362, 25.843, 2.490, 8.820 ], 
              [ 21.088, 1.357, 3.760, 26.029, 25.383, 1.926, 4.211, 25.843, 2.490, 20.000 ], 
              [ 22.596, 0.495, 22.782, 1.796, 2.500, 1.232, 1.800, 1.000, 1.000, 20.000 ], 
              [ 24.412, 1.567, 25.000, 26.020, 25.383, 25.000, 2.500, 25.843, 2.490, 20.000 ],
              [ 4.94, 0.15, 4.45, 15.96, 1.83, 1.96, 9.92, 9.33, 0.25, 39.32],
              [ 3.53, 0.54, 5.08, 16.61, 10.02, 6.25, 4.68, 1.03, 1.32, 12.14],
              [ 7.67, 0.51, 4.13, 15.57, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14], 
              [ 5.76, 0.54, 5.08, 9.64, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14],
              [ 4.94, 0.47, 4.77, 1.02, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14],
              [ 8.37, 0.54, 5.08, 9.64, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14],
              [ 7.14, 0.47, 3.10, 3.42, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14], 
              [ 5.22, 0.22, 4.45, 17.70, 2.97, 1.96, 9.92, 9.33, 0.25, 39.32],
              [ 3.56, 0.51, 4.86, 5.99, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14], 
              [ 6.77, 0.54, 5.08, 16.85, 10.02, 6.25, 4.68, 1.03, 1.32, 12.14],
              [ 8.09, 0.51, 3.65, 12.08, 10.02, 4.90, 8.88, 1.03, 1.32, 12.14], 
              [ 9.05, 0.49, 3.10, 1.02, 10.02, 4.90, 4.68, 1.03, 1.32, 12.14 ]]
# Synaptic factor values for selected Syns (0, 1, 3)
SynFactors =  [[4.5, 3.5, 20.0, 2.0, 0.001, 50.0, 5.8, 1.0, 0.5, 2.0, 2.5, 0.5], 
			  [ 4.0, 3.0, 15.0, 3.0, 0.0, 50.0, 1.5, 5.0, 2.0, 5.0, 2.0, 2.0], 
			  [ 3.5, 3.5, 7.0, 6.0, 0.04, 5.0, 6.5, 5.0, 4.5, 5.0, 8.0, 4.5], 
			  [ 1.5, 2.4, 4.0, 5.0, 0.04, 5.0, 4.5, 5.0, 4.5, 5.0, 5.0, 4.5]]
# Synaptic conductances for each synaptic input
SynConduct = [0.0002, 0.0002, 0.0002, 0.0002, 0.00115, 0.00051, 0.00052, 0.000041, 0.000065, 0.0003, 0.0002, 0.00037]
# Porportion of boutons that goes to somatoapical dendrite
PropBoutons = [[0.035, 0.000, 0.073, 0.097, 0.000, 0.000, 0.000, 0.013, 0.135, 0.154, 0.000, 0.000], 
              [ 0.305, 0.000, 0.595, 0.455, 0.000, 0.008, 0.362, 0.017, 0.077, 0.081, 0.333, 0.090], 
              [ 0.282, 0.000, 0.592, 0.450, 0.000, 0.000, 0.362, 0.006, 0.069, 0.074, 0.333, 0.071], 
              [ 0.259, 0.000, 0.584, 0.440, 0.000, 0.000, 0.362, 0.005, 0.053, 0.062, 0.333, 0.071]]
# Numpy array
IntrFactors = np.array(IntrFactors)
SynFactors = np.array(SynFactors)
SynConduct = np.array(SynConduct)
PropBoutons = np.array(PropBoutons)
# Index of CA3, EC3, InhSoma, InhDend
idxCA3 = 0
idxCA2 = 1
idxEC3 = 2
idxEC2 = 3
idxInhSoma = [4,6,10]
idxInhDend = [5,7,8,9,11]
# Index of Na, A, KDR, HCN, L, Ra
idxNa = 0
idxA = 1
idxAHPs = 2
idxC = 3
idxKDR = 4
idxM = 5
idxCa = 6
idxHCN = 7
idxL = 8
idxRa = 9
# Number of apical dendritic branches
NumApicDend = [ 48, 40, 44, 67 ]
# Number of basal dendritic branches
NumBasalDend = [ 46, 41, 22, 30 ]
# Morphology names
Morphologies = ['sup1','n409','n127','n128']




# 		==================================
# 			DEFINING INDIVIDUALS			 
# 		==================================

# Final 20 individuals
Individuals = [ [0,0], [0,1], [1,0], [1,1], [1,3],
				[2,0], [2,1], [4,0], [4,1], [6,1],
				[6,3], [7,1], [9,1], [12,3], [15,0],
				[15,3], [19,1], [19,3], [20,3], [21,0]]


# 		==================================
# 			DEFINING DATABASE	 
# 		==================================

Data = {  'Spikes': [[] for ii in range(len(Morphologies)*len(Individuals))],
		    'PrefPhase': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'Indiv': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'Intr': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'Syn': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'Morph': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'CA3': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'CA2': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'EC3': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'EC2': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'InhSoma': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'InhDend': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'ApicBranches': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'BasalBranches': [np.nan for ii in range(len(Morphologies)*len(Individuals))],
		    'TotalBranches': [np.nan for ii in range(len(Morphologies)*len(Individuals))] }



# 		==================================
# 			FILLING DATABASE	 
# 		==================================

numCycs = []
with open(DIR_LOCATION+'/LCNhm-results/TimeSpikes_by_individual.txt') as file:
	for line in file:
		Intr = int(line.split(' ')[0])
		Syn = int(line.split(' ')[1])
		Morph = line.split(' ')[2]
		# If it's a selected individual
		if (Morph in Morphologies) and ([Intr,Syn] in Individuals):
			iIndiv = Individuals.index([Intr,Syn])
			iMorph = Morphologies.index(Morph)
			iList = iIndiv*len(Morphologies) + iMorph
			
			# All Spikes
			if all(np.array(line[:-1].split(' ')[3:])==''): numCycs.append(len(line[:-1].split(' ')[3:]))
			Spks = np.array(map(float,[ii for ii in line[:-1].split(' ')[3:] if ii!='']))
			Spks = Spks[Spks>60.]
			Data['Spikes'][iList] = Spks
			
			# Preferred phase
			PrefPh = compute_prefer_phase(Spks)
			Data['PrefPhase'][iList] = PrefPh

			# Number of Intrinsic and Syaptic set, and Morphology
			Data['Indiv'][iList] = iIndiv
			Data['Intr'][iList] = Intr
			Data['Syn'][iList] = Syn
			Data['Morph'][iList] = Morph

			# Absolute conductance of CA3 and EC3, as: 
			# 	(% Boutons in somatoapical trunk) x (Synaptic Factor) x (Synaptic Conductance)
			Data['CA3'][iList] = PropBoutons[iMorph,idxCA3] * SynFactors[Syn,idxCA3] * SynConduct[idxCA3]
			Data['EC3'][iList] = PropBoutons[iMorph,idxEC3] * SynFactors[Syn,idxEC3] * SynConduct[idxEC3]
			Data['CA2'][iList] = 1.0 * SynFactors[Syn,idxCA2] * SynConduct[idxCA2]
			Data['EC2'][iList] = 1.0 * SynFactors[Syn,idxEC2] * SynConduct[idxEC2]
			Data['InhSoma'][iList] = np.sum([PropBoutons[iMorph,idxInh] * SynFactors[Syn,idxInh] * SynConduct[idxInh] for idxInh in idxInhSoma])
			Data['InhDend'][iList] = np.sum([PropBoutons[iMorph,idxInh] * SynFactors[Syn,idxInh] * SynConduct[idxInh] for idxInh in idxInhDend])
			Data['ApicBranches'][iList] = NumApicDend[iMorph]
			Data['BasalBranches'][iList] = NumBasalDend[iMorph]
			Data['TotalBranches'][iList] = NumApicDend[iMorph] + NumBasalDend[iMorph]



# 		==================================
# 			PLOT ABSOLUTE SYNAPTIC	 
# 		==================================

ys = np.array(Data['PrefPhase'])
xNames = ['CA3','CA2','EC3','EC2','InhSoma','InhDend']
plt.figure(figsize=(len(xNames)*3,3))
for ii, var in enumerate(xNames):
	plt.subplot(1,len(xNames),ii+1)
	xs = np.array(Data[var])
	xsUniq = np.unique(xs)
	ysUp = []
	ysDown = []
	for xUniq in xsUniq:
		ysUp.append(np.mean(ys[xs==xUniq])+np.std(ys[xs==xUniq]))
		ysDown.append(np.mean(ys[xs==xUniq])-np.std(ys[xs==xUniq]))

	plt.fill_between(xsUniq, ysDown, ysUp, color='g', alpha=0.3)
	plt.scatter(xs, ys, color='k', alpha=0.2, s=40)
	pval, c, a, R2 = compute_pval(xs, ys)
	plt.plot(xsUniq, c + a*xsUniq, '--g')
	plt.title('R = %.2f\npval = %.4f'%(np.sqrt(R2), pval))
	plt.xlabel(var)
	plt.yticks(np.arange(0,361,90))
	plt.xlim([min(xs)-0.1*max(xs),1.1*max(xs)])
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_SynapticAbsolute.png', bbox_inches="tight")
#plt.show()



# 		==================================
# 			PLOT SYNAPTIC CONDUCTANCE
# 		==================================

ys = np.array(Data['PrefPhase'])
xNames = ['CA3','CA2','EC3','EC2','InhSoma','InhDend']
xVars = [idxCA3,idxCA2,idxEC3,idxEC2,idxInhSoma,idxInhDend]
plt.figure(figsize=(len(xNames)*3,3))
for ii, var in enumerate(xNames):
	plt.subplot(1,len(xNames),ii+1)
	if type(xVars[ii])==int:
		xs = np.array(SynFactors[Data['Syn'],xVars[ii]])
	else:
		xs = np.sum([SynFactors[Data['Syn'],jj] for jj in xVars[ii]],0)
	xsUniq = np.unique(xs)
	ysUp = []
	ysDown = []
	for xUniq in xsUniq:
		ysUp.append(np.mean(ys[xs==xUniq])+np.std(ys[xs==xUniq]))
		ysDown.append(np.mean(ys[xs==xUniq])-np.std(ys[xs==xUniq]))

	plt.fill_between(xsUniq, ysDown, ysUp, color='g', alpha=0.3)
	plt.scatter(xs, ys, color='k', alpha=0.2, s=40)
	pval, c, a, R2 = compute_pval(xs, ys)
	plt.plot(xsUniq, c + a*xsUniq, '--g')
	plt.title('R = %.2f\npval = %.4f'%(np.sqrt(R2), pval))
	plt.xlabel(var)
	plt.yticks(np.arange(0,361,90))
	plt.xlim([min(xs)-0.1*max(xs),1.1*max(xs)])
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_SynapticConductance.png', bbox_inches="tight")
#plt.show()


# 		==================================
# 			PLOT INTRINSIC	 
# 		==================================

ys = np.array(Data['PrefPhase'])
xNames = ['idxNa','idxA','idxAHPs','idxC','idxKDR','idxM','idxCa','idxHCN','idxL','idxRa']
xVars = [idxNa,idxA,idxAHPs,idxC,idxKDR,idxM,idxCa,idxHCN,idxL,idxRa]
plt.figure(figsize=(len(xNames)*3,3))
for ii, var in enumerate(xNames):
	plt.subplot(1,len(xNames),ii+1)
	xs = np.array(IntrFactors[Data['Intr'],xVars[ii]])
	xsUniq = np.unique(xs)
	ysUp = []
	ysDown = []
	for xUniq in xsUniq:
		ysUp.append(np.mean(ys[xs==xUniq])+np.std(ys[xs==xUniq]))
		ysDown.append(np.mean(ys[xs==xUniq])-np.std(ys[xs==xUniq]))

	pval, c, a, R2 = compute_pval(xs, ys)
	plt.fill_between(xsUniq, ysDown, ysUp, color='g', alpha=0.3)
	plt.scatter(xs, ys, color='k', alpha=0.2, s=40)
	plt.plot(xsUniq, c + a*xsUniq, '--g')
	plt.title('R = %.2f\npval = %.4f'%(np.sqrt(R2), pval))
	plt.xlabel(var[3:])
	plt.yticks(np.arange(0,361,90))
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_Intrinsic.png', bbox_inches="tight")
#plt.show()



# 		==================================
# 			PLOT MORPHOLOGY	 
# 		==================================

ys = np.array(Data['PrefPhase'])
xNames = ['ApicBranches','BasalBranches','TotalBranches']
plt.figure(figsize=(len(xNames)*3,3))
for ii, var in enumerate(xNames):
	plt.subplot(1,len(xNames),ii+1)
	xs = np.array(Data[var])
	xsUniq = np.unique(xs)
	ysUp = []
	ysDown = []
	for xUniq in xsUniq:
		ysUp.append(np.mean(ys[xs==xUniq])+np.std(ys[xs==xUniq]))
		ysDown.append(np.mean(ys[xs==xUniq])-np.std(ys[xs==xUniq]))

	pval, c, a, R2 = compute_pval(xs, ys)
	plt.fill_between(xsUniq, ysDown, ysUp, color='g', alpha=0.3)
	plt.scatter(xs, ys, color='k', alpha=0.2, s=40)
	plt.plot(xsUniq, c + a*xsUniq, '--g')
	plt.title('R = %.2f\npval = %.4f'%(np.sqrt(R2), pval))
	plt.xlabel(var)
	plt.yticks(np.arange(0,361,90))
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_Morphology.png', bbox_inches="tight")
plt.show()



# 		==================================
# 			PLOT INTR+SYN BOXPLOT
# 		==================================


# Matrix of values for each individual
intrNames = ['gNa','gA','gKDR','gHCN','gL']
intrVars = [idxNa,idxA,idxKDR,idxHCN,idxL]
intrMat = np.zeros((len(Individuals),len(intrNames)))
for iVar, intrName in enumerate(intrNames):
	intrMat[:,iVar] = IntrFactors[[ii for ii,jj in Individuals],varVars[iVar]]
intrMat = np.array(intrMat)
intrMatNorm = (intrMat-np.min(intrMat,0))/(np.max(intrMat,0)-np.min(intrMat,0))
# Matrix of values for each individual
synNames = ['gCA3','gCA2','gEC3','gInhSoma','gInhDend']
synVars = [idxCA3,idxCA2,idxEC3,idxInhSoma,idxInhDend]
synMat = np.zeros((len(Individuals),len(synNames)))
for iVar, synName in enumerate(synNames):
	if type(synVars[iVar])==int:
		synMat[:,iVar] = SynConduct[synVars[iVar]] * SynFactors[[jj for ii,jj in Individuals],synVars[iVar]]
	else:
		synMat[:,iVar] = np.sum([SynConduct[kk] * SynFactors[[jj for ii,jj in Individuals],kk] for kk in synVars[iVar]],0)
synMat = np.array(synMat)
synMatNorm = (synMat-np.min(synMat,0))/(np.max(synMat,0)-np.min(synMat,0))

# Prefered phases for each individual
PrefPhases = np.array(Data['PrefPhase'])
PrefPhases[PrefPhases<50] = PrefPhases[PrefPhases<50]+360
Indivs = np.array(Data['Indiv'])
PrefPhIndivs = np.array([ PrefPhases[Indivs==iI] for iI in range(len(Individuals))])
PrefPhIndivMean = [ np.mean(PrefPhases[Indivs==iI]) for iI in range(len(Individuals))]
PrefPhIndivMedian = [ np.median(PrefPhases[Indivs==iI]) for iI in range(len(Individuals))]
PrefPhIndivStd = [ np.std(PrefPhases[Indivs==iI]) for iI in range(len(Individuals))]
PrefPhOrder = np.argsort(-np.array(PrefPhIndivMedian))

# Figure
plt.figure(figsize=(12,12))
gs = gridspec.GridSpec(1, 3, width_ratios=[3, 1, 1]) 
# Phases Boxplot 
ax0 = plt.subplot(gs[0])
ax0.boxplot(np.transpose(PrefPhIndivs[PrefPhOrder]),vert=False,showmeans=True,meanline=True, showfliers=False, notch=True)
Hight = 0
for iIndiv in PrefPhOrder:
	ax0.scatter(PrefPhases[Indivs==iIndiv], [Hight+1]*4, 10, color='k')
	#ax0.text(-15,Hight+1,iIndiv)
	ax0.text(-15,Hight+1,'(%d,%d)'%tuple(Individuals[iIndiv]))
	Hight += 1
ax0.plot( np.arange(0,450), 19/3.*(1-np.cos(np.arange(0,450.)/180.*np.pi)),'k',linewidth=2)
ax0.set_yticks([])
ax0.set_xticks(np.arange(0,451,180))
ax0.set_ylim([-2,len(Individuals)+2])
# Intrinsic Values
ax1 = plt.subplot(gs[1])
ax1.imshow(np.exp(3*(1-intrMatNorm[PrefPhOrder])),cmap='summer')
ax1.set_yticks([])
ax1.set_xticks(range(len(intrNames)))
ax1.set_xticklabels(intrNames,rotation=-90)
# Synaptic Values
ax2 = plt.subplot(gs[2])
ax2.imshow(np.exp(3*(1-synMatNorm[PrefPhOrder])),cmap='summer')
ax2.set_yticks([])
ax2.set_xticks(range(len(synNames)))
ax2.set_xticklabels(synNames,rotation=-90)

plt.tight_layout()
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_IntrSynMat.png', bbox_inches="tight")
plt.show()



# 		==================================
# 			PLOT MORPHOLOGY BOXPLOT
# 		==================================


# Matrix of values for each morphology
morphNames = ['ApicBranches','BasalBranches','TotalBranches','PropBoutCA3','PropBoutEC3','PropBoutInhSoma','PropBoutInhDend']
morphVars = [NumApicDend,NumBasalDend,np.sum([NumApicDend,NumBasalDend],0),PropBoutons[:,idxCA3],PropBoutons[:,idxEC3],np.sum(PropBoutons[:,idxInhSoma],1),np.sum(PropBoutons[:,idxInhDend],1)]
morphMat = np.zeros((len(Morphologies),len(morphNames)))
for iVar, morphName in enumerate(morphNames):
	morphMat[:,iVar] = morphVars[iVar]
morphMat = np.array(morphMat)
morphMatNorm = (morphMat-np.min(morphMat,0))/(np.max(morphMat,0)-np.min(morphMat,0))

# Prefered phases for each individual
PrefPhases = np.array(Data['PrefPhase'])
PrefPhases[PrefPhases<50] = PrefPhases[PrefPhases<50]+360
Morphs = np.array(Data['Morph'])
PrefPhMorphs = np.array([ PrefPhases[Morphs==Morph] for Morph in Morphologies])
PrefPhMorphMean = [ np.mean(PrefPhases[Morphs==Morph]) for Morph in Morphologies]
PrefPhMorphMedian = [ np.median(PrefPhases[Morphs==Morph]) for Morph in Morphologies]
PrefPhMorphStd = [ np.std(PrefPhases[Morphs==Morph]) for Morph in Morphologies]
PrefPhOrder = np.argsort(-np.array(PrefPhMorphMedian))

# Figure
plt.figure(figsize=(12,12))
gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1]) 
# Phases Boxplot 
ax0 = plt.subplot(gs[0])
ax0.boxplot(np.transpose(PrefPhMorphs[PrefPhOrder]),vert=False, showmeans=True, meanline=True, showfliers=False, notch=True)
Hight = 0
for iMorph in PrefPhOrder:
	ax0.scatter(PrefPhases[Morphs==Morphologies[iMorph]], [Hight+0.75]+np.random.rand(20)*0.5, 10, color='k')
	ax0.text(-15,Hight+1,Morphologies[iMorph])
	Hight += 1
ax0.plot( np.arange(0,450), 4/3.*(1-np.cos(np.arange(0,450.)/180.*np.pi)),'k',linewidth=2)
ax0.set_yticks([])
ax0.set_xticks(np.arange(0,451,180))
ax0.set_ylim([-0.1,len(Morphologies)+.6])
# Intrinsic Values
ax1 = plt.subplot(gs[1])
ax1.imshow(np.exp(2.2*(1-morphMatNorm[PrefPhOrder])),cmap='summer',aspect='auto')
ax1.set_yticks([])
ax1.set_xticks(range(len(morphNames)))
ax1.set_xticklabels(morphNames,rotation=-90)

plt.tight_layout()
plt.savefig(DIR_LOCATION+'/LCNhm-results/thetaBasal_statistics_MorphMat.png', bbox_inches="tight")
plt.show()

