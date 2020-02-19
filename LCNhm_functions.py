"""
This script has all functions except from the :ref:`LCNhm-class` and its methods.

Within its miscellaneous purposes, there are:

* One function naming and making the folder: ``make_folder``

* Functions related to how to distribute intrinsic properties along the cell: ``gradient_membrane_resistance``, ``gradient_spine_scale``, ``gradient_V12`` and ``gradient_gHCN``

* Function of how to distribute the synaptic releases along a theta cycle: ``synaptic_time_probability_distribution``

* Functions that record and save data from the simulation: ``recordings``, ``save_spiking_times``, ``save_recordings`` and ``save_parameters``

Functions and descriptions
--------------------------
"""
import os
import sys
import datetime
import numpy as np
from neuron import h, nrn, gui

def make_folder(DIR_LOCATION, OPT_FOLDER_NAME):
    """
    Make a new folder inside *DIR_LOCATION/LCNhm-results/* to save any new simulation.

        ``< D(ate)YearMonthDay _ T(ime)HourMinute _ OptionalFolderName >``

    Parameters
    ----------
    DIR_LOCATION: String
        Name of main directory, where *LCNhm_main.py* is

        Defined in :const:`LCNhm_configurationfile.DIR_LOCATION`

    OPT_FOLDER_NAME: String
        Optional folder name: additional description of the simulation.

        Defined in :const:`LCNhm_configurationfile.OPT_FOLDER_NAME`
        
        Eg.: *first_test*, *high_CA3*, *inhibition_suppression*

    Returns
    -------
    FOLDER_NAME: String
        Name of the folder just made.
    """

    # Today's date and time
    Date = [ datetime.datetime.today().year,
             datetime.datetime.today().month,
             datetime.datetime.today().day,
             datetime.datetime.today().hour,
             datetime.datetime.today().minute ]

    # Name of folder: <D(ate)YearMonthDay_T(ime)HourMinute_OptionalFolderName>
    FOLDER_NAME = 'D%.2d%.2d%.2d_T%.2d%.2d'%tuple(Date) + '_' * (len(OPT_FOLDER_NAME)>0) + OPT_FOLDER_NAME
    FolderLen = len(FOLDER_NAME)

    # Make a new directory
    dirFolderName = DIR_LOCATION+'/LCNhm-results/'+FOLDER_NAME

    # If folder doesnt exist, make a new one; if it exists, add a suffix
    k = 1
    while os.path.exists(dirFolderName):
        FOLDER_NAME = FOLDER_NAME[:FolderLen] + '_%.3d'%k
        dirFolderName = DIR_LOCATION+'/LCNhm-results/'+FOLDER_NAME
        k += 1
    # Make new folder
    os.makedirs(dirFolderName)

    return dirFolderName

def gradient_membrane_resistance(zDistance):
    """
    Definition of the membrane resistivity along the *z* axis (radial axis)

    Parameters
    ----------
    zDistance: Float
        Distance to soma on the *z* axis (radial axis)

    Returns
    -------
    MembraneResistance: Float
        Membrane resistance (ohm*cm2) at :data:`zDistance`

    """
    # Conductance:  gm = 1/Rm   [S/cm2]
    #               Rm = 1/gm   [Ohm*cm2]
    kohm2ohm = 1000
    MembraneResistance = kohm2ohm*(80. + (0.4-80.)/(1+np.exp((225-zDistance)/30.)))   # ohm * cm2
    return MembraneResistance

def gradient_spine_scale(zDistance, Diameter):
    """
    Definition of the spine scale factor along the *z* axis (radial axis)

    To take into account the spines intrinsic properties such as the 
    capacitance and the membrane resistance must be multiplied by a factor SS

    Parameters
    ----------
    zDistance: Float
        Distance to soma on the *z* axis (radial axis)

    Diameter: Float
        Diameter of compartiment (um)

    Returns
    -------
    SS: Float
        Spine scale factor            
    """
    
    # Str Lacunosum Moleculare
    if zDistance>350:
        if (Diameter>0.35): SS = 2.10 # Medium
        else: SS = 1.71 # Thin
    # Str Radiatum
    elif (zDistance<= 350) and (zDistance>0):
        if Diameter>1.6: SS = 1.69 # Thick medial
        elif (Diameter<=1.6) and (Diameter>0.55): SS = 1.60 # Thick distal
        elif (Diameter<0.55): SS = 1.86 # Thin
        else: SS = 1
    # Str Oriens
    elif zDistance<0: SS = 2.51
    else: SS = 1

    return SS

def gradient_V12(zDistance):
    """
    Definition of the HCN half-maximal activation voltage (V1/2) factor along the *z* axis (radial axis) (Sinha2015)

    Parameters
    ----------
    zDistance: Float
        Distance to soma on the *z* axis (radial axis)

    Returns
    -------
    v12: Float
        Half-maximal activation voltage (mV)
    """
    v12 = -82*(zDistance<300) - (0.04*zDistance-0.04*100)*(zDistance>100 and zDistance<300) - 90*(zDistance>=300)   # mV
    return v12

def gradient_gHCN(zDistance):
    """
    Definition of the HCN maximum conductance along the *z* axis (radial axis) (Sinha2015)

    Parameters
    ----------
    zDistance: Float
        Distance to soma on the *z* axis (radial axis)

    Returns
    -------
    gHCN: Float
        Maximumn HCN conductance (S/cm2)
    """
    ghcnBase = 85.e-6   # S/cm2
    gHCN = ghcnBase*(1+120./(1. + np.exp((275.-zDistance)/50.)) )   # S/cm2
    return gHCN

def synaptic_time_probability_distribution(Times, PhaseMax, Parameters):
    """
    Definition of the synaptic probability distribution along theta phase.
    In this case, the distribution chosen is the asymmetric gaussian (beta function),
    described as:

    .. code-block:: python

        beta(x,A,B) = x^(A-1) * (1-x)^(B-1)

    Parameters
    ----------
    Times: Numpy array
        Times to whom compute the firing probability function

    PhaseMax: Float
        Phase of maximum probability

    Parameters: List
        List of parameters. 
        In this case, ``Parameters = [A, B]`` will define the shift of the gaussian distribution
        to the left or to the right.

    Returns
    -------
    ProbDistribution: Numpy array
        Spiking probability distribution for ``Times``
    """

    # Unpack parameters
    ThetaPeriod, A, B = Parameters

    # Phases
    Phases = np.mod(Times,ThetaPeriod)/ThetaPeriod
    dPh = Phases[1]-Phases[0]
    # Phase (desired position of the maximum)
    PhaseMax = np.argwhere(np.abs(Phases-PhaseMax/360.)<=(Phases[1]-Phases[0]))[0,0]
    # Actual position of the maximum
    Phasesmax = np.argmax(Phases**(A-1.)*(1-Phases)**(B-1.))
    while dPh*(Phasesmax-PhaseMax)>1:
        Phasesmax -= int(1./dPh)

    # Shifting the distribution
    idx = 1+dPh*(Phasesmax-PhaseMax)
    Phases2 = np.mod(np.array([Phases[-1]+i*dPh for i in np.arange(int(idx/dPh))]) , 1)
    Phases = Phases[int(idx/dPh):]

    ProbDistribution = np.append(Phases**(A-1.)*(1-Phases)**(B-1.), Phases2**(A-1.)*(1-Phases2)**(B-1.) ) 
    ProbDistribution /= sum(ProbDistribution)    

    return ProbDistribution

def recordings(Pyramidal, RECORDING_MAGNITUDE, RECORDING_SECTION, RECORDING_LOCATION):
    """
    Recordings are set.

    Parameters
    ----------
    Pyramidal: :class:`LCNhm_class.neuron_class` object

    RECORDING_MAGNITUDE: List of `str`
        Physical magnitudes to be measured and recorded

        Defined in :const:`LCNhm_configurationfile.RECORDING_MAGNITUDE`

    RECORDING_SECTION: List of `str`
        Sites/Sections to be measured and recorded

        Defined in :const:`LCNhm_configurationfile.RECORDING_SECTION`

    RECORDING_LOCATION: List of Floats
        Location along the defined :const:`RECORDING_SECTION` to be measured and recorded

        Defined in :const:`LCNhm_configurationfile.RECORDING_LOCATION`

    Returns
    -------
    Recordings: Dictionary

        Dictionary with all recordings

            Eg.: Given the following inputs, the Recordings Dictionary would be
            
            .. code-block:: python

                RECORDING_MAGNITUDE = ['Time','Vmem','Pos']
                RECORDING_SECTION = ['SomaList0','ApicList0','ApicList0','ApicList0']
                RECORDING_LOCATION = [0.0, 0.2, 0.5, 0.9]
                
                Recordings = { 'Time': record object
                               'Vmem': { 'SomaList0_000': record object, 
                                         'ApicList0_020': record object, 
                                         'ApicList0_050': record object, 
                                         'ApicList0_090': record object}, 
                               'Imem': { 'SomaList0_000': record object, 
                                         'ApicList0_020': record object, 
                                         'ApicList0_050': record object, 
                                         'ApicList0_090': record object} }

    """
    # If it's empty, record Soma for Spike times, but don't save it
    if len(RECORDING_MAGNITUDE) == 0:
        RECORDING_MAGNITUDE = ['Vmem']
        RECORDING_SECTION = ['SomaList0']
        RECORDING_LOCATION = [0.0]

    # Initialize variable
    Recordings = {}

    for Mag in RECORDING_MAGNITUDE:
        Recordings[Mag] = {}

    # Time
    Recordings['Time'] = h.Vector()
    Recordings['Time'].record( h._ref_t )
    # Magnitudes to be recorded
    for Mag in RECORDING_MAGNITUDE:
        # Sections defined in RECORDING_SECTION
        for ii, Section in enumerate(RECORDING_SECTION):
            # Get section object
            sec = getattr(Pyramidal,Section[:8])[int(Section[8:])]
            loc = RECORDING_LOCATION[ii]
            seg = np.floor(loc*sec.nseg)
            SecName = '%s_%.3d'%(Section,loc*100)
            if Mag!='Time': 
                Recordings[Mag][SecName] = h.Vector()
                if Mag=='Vmem': Recordings[Mag][SecName].record( sec( loc )._ref_v )
                if Mag=='Imem': Recordings[Mag][SecName].record( sec( loc )._ref_i_membrane )
                if Mag=='Pos': Recordings[Mag][SecName] = [h.x3d(seg,sec=sec), h.y3d(seg,sec=sec), h.z3d(seg,sec=sec), h.diam3d(seg,sec=sec)]

    # Always record Vmem of SomaList0 and Time:
    Mag = 'Vmem'
    SecName = 'SomaList0_000'
    if Mag not in Recordings.keys():
        Recordings[Mag] = {}
        Recordings[Mag][SecName] = h.Vector()
        Recordings[Mag][SecName].record( sec( loc )._ref_v )
    elif SecName not in Recordings[Mag].keys():
        Recordings[Mag][SecName] = h.Vector()
        Recordings[Mag][SecName].record( sec( loc )._ref_v )

    return Recordings

def save_spiking_times(Recordings, FolderName):
    """
    Time of spikes are detected and saved in :data:`LCNhm_main.FolderName` */TimeSpikes.txt*

    If there are no spikes, it will write ``NaN``

    Parameters
    ----------
    Recordings: Dictionary
        Dictionary with all recordings, output from :func:`recordings`

    FolderName: String
        Name of folder where recordings will be saved, output from :func:`make_folder`
    """

    # Membrane potential of soma and time
    VmemSoma = np.array(Recordings['Vmem']['SomaList0_000'].to_python())
    Time = np.array(Recordings['Time'].to_python())

    # Initialize variable
    TimeSpikes = []
    # Time of spikes
    for i, iTime in enumerate(Time[Time>50.]):
        if (VmemSoma[i]>VmemSoma[i-1]) and (VmemSoma[i]>VmemSoma[i+1]) and (VmemSoma[i]>-10):
            TimeSpikes.append(iTime)

    # If there are no spikes, 'nan'
    if len(TimeSpikes)==0:
        TimeSpikes.append(np.nan)

    # Save
    #with open('/'.join(FolderName.split('/')[:-1])+'/TimeSpikes.txt',"a") as file:
    with open(FolderName+'/TimeSpikes.txt',"a") as file:
        file.write(" ".join(map(str, TimeSpikes ))) 
        file.write("\n")

def save_recordings(Recordings, FolderName, RECORDING_MAGNITUDE):
    """
    Recordings of the magnitudes selected in :const:`LCNhm_configurationfile.RECORDING_MAGNITUDE`
    are saved in a txt file in :data:`LCNhm_main.FolderName` subfolder, as:

        Eg.: *Recordings_Vmem.txt*

        Eg.: *Recordings_Pos.txt*

    Membrane potentials (Vmem) and membrane currents (Imem) along time of each recording site 
    is saved in different rows

        Eg.:

        >>> vi LCNhm-results/20190723_1200_test/Recordings_Vmem.txt
        >>>           Time    0.000    0.025    0.050   ...  
        >>>  SomaList0_000  -65.000  -65.300  -65.450   ...  
        >>>  ApicList0_000  -55.000  -54.800  -54.600   ...  

    Position is saved in different rows for each recording site 

        Eg.:
        
        >>> vi LCNhm-results/20190723_1200_test/Recordings_Pos.txt
        >>>  SomaList0_000  0.0    0.0      0.0
        >>>  ApicList0_000  0.0    0.0  -1000.0


    Parameters
    ----------
    Recordings: Dictionary
        Dictionary with all recordings, output from :func:`recordings`

    FolderName: String
        Name of folder where recordings will be saved, output from :func:`make_folder`
    """
    # If RECORDING_SECTION is empty, record Soma for Spike times, but don't save it
    if len(RECORDING_MAGNITUDE)>0:
        # Magnitudes to be recorded
        for Mag in Recordings.keys():
            if Mag!='Time':
                if Mag=='Pos':
                    with open('%s/Recordings_%s.txt'%(FolderName,Mag),"a") as file:
                        # Sections
                        for SecName in Recordings[Mag].keys():
                            # Write magnitude of section 'SecName'
                            file.write("%s %.3f %.3f %.3f %.3f\n"%tuple([SecName]+Recordings[Mag][SecName])) 
                else:
                    with open('%s/Recordings_%s.txt'%(FolderName,Mag),"a") as file:
                        # Time
                        file.write("Time "+" ".join(map(str, Recordings['Time'].to_python() ))) 
                        file.write("\n")
                        # Sections
                        for SecName in Recordings[Mag].keys():
                            # Write magnitude of section 'SecName'
                            file.write(SecName+" "+" ".join(map(str, Recordings[Mag][SecName].to_python() ))) 
                            file.write("\n")

def save_parameters(Parameters, FolderName):
    """
    Values of the :ref:`LCNhm-configuration-file` parameters are written in Parameters.txt
    in :data:`LCNhm_main.FolderName` subfolder.

        Eg.:

        >>> vi LCNhm-results/20190723_1200_test/Parameters.txt
        >>>   DIR_LOCATION /home/andrea/Projects/HippoModel/LCNhippomodel
        >>>   OPT_FOLDER_NAME test
        >>>   SIMPROP_THETA_MODE True
        >>>   SIMPROP_THETA_PERIOD 166.000000
        >>>   SIMPROP_START_TIME 0.000000
        >>>   SIMPROP_SIM_TIME 1660.000000
        >>>   SIMPROP_END_TIME 1660.000000
        >>>   SIMPROP_DT 0.025000
        >>>   SIMPROP_TEMPERATURE 34.000000
        >>>   CELLPROP_MORPHOLOGY n128
        >>>   CELLPROP_INTRINSIC 0
        >>>   CELLPROP_INTRINSIC_IONCHS iNas iA iAHPs iC iCaL iCaT iKDR iM iHCN iL
        >>>   CELLPROP_INTRINSIC_EXPERIMENT 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
        >>>   CELLPROP_SYNAPTIC_INPUTS CA3 CA2 EC3 EC2 Axo Bis CCK Ivy NGF OLM PV SCA
        >>>   CELLPROP_SYNAPTIC_EXPERIMENT 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000 1.000000
        >>>   CURRENT_DURATION 0.000000
        >>>   CURRENT_DELAY 0.000000
        >>>   CURRENT_AMPLITUDES 0.000000
        >>>   CURRENT_SECTION SomaList0
        >>>   CURRENT_LOCATION 0.000000
        >>>   RECORDING_MAGNITUDE 
        >>>   RECORDING_SECTION SomaList0
        >>>   RECORDING_LOCATION 0.0 

    Parameters
    ----------
    Parameters: List
        List with all :ref:`LCNhm-configuration-file` parameters

    FolderName: String
        Name of folder where recordings will be saved, output from :func:`make_folder`
    """

    # Save configuration file
    DIR_LOCATION, OPT_FOLDER_NAME, SIMPROP_THETA_MODE, SIMPROP_THETA_PERIOD, SIMPROP_START_TIME, SIMPROP_SIM_TIME, SIMPROP_END_TIME, SIMPROP_DT, SIMPROP_TEMPERATURE, CELLPROP_MORPHOLOGY, CELLPROP_INTRINSIC, CELLPROP_SYNAPTIC, CELLPROP_INTRINSIC_IONCHS, CELLPROP_INTRINSIC_EXPERIMENT, CELLPROP_SYNAPTIC_INPUTS, CELLPROP_SYNAPTIC_EXPERIMENT, CURRENT_DURATION, CURRENT_DELAY, CURRENT_AMPLITUDES, CURRENT_SECTION, CURRENT_LOCATION, RECORDING_MAGNITUDE, RECORDING_SECTION, RECORDING_LOCATION = Parameters

    with open('%s/Parameters.txt'%FolderName,"a") as file:
        file.write('DIR_LOCATION %s\n'%DIR_LOCATION)
        file.write('OPT_FOLDER_NAME %s\n'%OPT_FOLDER_NAME)
        file.write('SIMPROP_THETA_MODE %s\n'%SIMPROP_THETA_MODE)
        file.write('SIMPROP_THETA_PERIOD %f\n'%SIMPROP_THETA_PERIOD)
        file.write('SIMPROP_START_TIME %f\n'%SIMPROP_START_TIME)
        file.write('SIMPROP_SIM_TIME %f\n'%SIMPROP_SIM_TIME)
        file.write('SIMPROP_END_TIME %f\n'%SIMPROP_END_TIME)
        file.write('SIMPROP_DT %f\n'%SIMPROP_DT)
        file.write('SIMPROP_TEMPERATURE %f\n'%SIMPROP_TEMPERATURE)
        file.write('CELLPROP_MORPHOLOGY %s\n'%CELLPROP_MORPHOLOGY)
        file.write('CELLPROP_INTRINSIC %d\n'%CELLPROP_INTRINSIC)
        file.write('CELLPROP_SYNAPTIC %d\n'%CELLPROP_SYNAPTIC)
        file.write('CELLPROP_INTRINSIC_IONCHS'+' %s'*len(CELLPROP_INTRINSIC_IONCHS)%tuple(CELLPROP_INTRINSIC_IONCHS)+'\n')
        file.write('CELLPROP_INTRINSIC_EXPERIMENT'+' %f'*len(CELLPROP_INTRINSIC_EXPERIMENT)%tuple(CELLPROP_INTRINSIC_EXPERIMENT)+'\n')
        file.write('CELLPROP_SYNAPTIC_INPUTS'+' %s'*len(CELLPROP_SYNAPTIC_INPUTS)%tuple(CELLPROP_SYNAPTIC_INPUTS)+'\n')
        file.write('CELLPROP_SYNAPTIC_EXPERIMENT'+' %f'*len(CELLPROP_SYNAPTIC_EXPERIMENT)%tuple(CELLPROP_SYNAPTIC_EXPERIMENT)+'\n')
        file.write('CURRENT_DURATION'+' %f'*len(CURRENT_DURATION)%tuple(CURRENT_DURATION)+'\n')
        file.write('CURRENT_DELAY'+' %f'*len(CURRENT_DELAY)%tuple(CURRENT_DELAY)+'\n')
        file.write('CURRENT_AMPLITUDES'+' %f'*len(CURRENT_AMPLITUDES)%tuple(CURRENT_AMPLITUDES)+'\n')
        file.write('CURRENT_SECTION'+' %s'*len(CURRENT_SECTION)%tuple(CURRENT_SECTION)+'\n')
        file.write('CURRENT_LOCATION'+' %f'*len(CURRENT_LOCATION)%tuple(CURRENT_LOCATION)+'\n')
        file.write('RECORDING_MAGNITUDE \n'%RECORDING_MAGNITUDE)
        file.write('RECORDING_SECTION'+' %s'*len(RECORDING_SECTION)%tuple(RECORDING_SECTION)+'\n')
        file.write('RECORDING_LOCATION'+' %s'*len(RECORDING_LOCATION)%tuple(RECORDING_LOCATION)+'\n')

