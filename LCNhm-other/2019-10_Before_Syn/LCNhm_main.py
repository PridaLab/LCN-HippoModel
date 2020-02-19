"""
It is the main script, the core module of the LCN-HippoModel project.

Its function is to:

1. Load the parameters from the :ref:`LCNhm-configuration-file`

#. Create the pyramidal cell according to the configuration parameters through the :ref:`LCNhm-class`

#. Run the simulation

#. Save the specified physical magnitudes inside ``DIR_LOCATION``/LCNhm-results/``FolderName``

The ``FolderName`` is define by one of the :ref:`LCNhm-functions`, and looks like:

	``< D(ate)YearMonthDay _ T(ime)HourMinute _ OPT_FOLDER_NAME >``

	Eg.: ./ LCNhm-results / D20190730_T1030_First_Test

It will add suffixes if the folder already exists:

	Eg.: ./ LCNhm-results / D20190730_T1030_First_Test_001
	
	Eg.: ./ LCNhm-results / D20190730_T1030_First_Test_002
	
	Eg.: ./ LCNhm-results / D20190730_T1030_First_Test_003
	
	Eg.: ./ LCNhm-results / D20190730_T1031_First_Test
	
	Eg.: ./ LCNhm-results / D20190730_T1031_First_Test_001

It will contain:

* **TimeSpikes.txt**: Time of spikes. If there are no spikes, it will write ``NaN``

* **Recordings_Vmem.txt**: Recordings of the membrane potential (if specified in the :ref:`LCNhm-configuration-file`)

* **Recordings_Imem.txt**: Recordings of the membrane current (if specified in the :ref:`LCNhm-configuration-file`)

* **Recordings_Pos.txt**: Spatial positions (if specified in the :ref:`LCNhm-configuration-file`)

Membrane potentials (Vmem) and membrane currents (Imem) along time of each recording site
is saved in different rows. Position is saved in different rows for each recording SIMPROP_TEMPERATURE


Functions and descriptions
--------------------------
"""
print '''
	------------------
	| LCN-HippoModel |
	------------------
'''

import sys
import time
import numpy as np
import pandas as pd
from neuron import h, nrn, gui

from LCNhm_configurationfile import *
from LCNhm_functions import *
from LCNhm_class import *

TimeCountStart = time.time()
# 	-----------
# 	Make folder
# 	-----------
FolderName = make_folder( DIR_LOCATION, OPT_FOLDER_NAME)
"""
String:
Full name of folder where results will be saved
"""
print '  %3d:%.2d ... folder "%s" made'%((time.time()-TimeCountStart)/60.,round((time.time()-TimeCountStart)%60),FolderName.split('/')[-1])

save_parameters(FolderName)

# 	---------------------------
# 	Make cell from neuron_class
# 	---------------------------
# Make inputs to class
CurrentFactors = [ CURRENT_DURATION, CURRENT_DELAY, CURRENT_AMPLITUDES,
				   CURRENT_SECTION, CURRENT_LOCATION ]
"""
List of properties for all desired current pulses, packed to be a one of the :class:`LCNhm_class.neuron_class`'s inputs

* :const:`LCNhm_configurationfile.CURRENT_DURATION` : List of durations of each current pulse in milliseconds

* :const:`LCNhm_configurationfile.CURRENT_DELAY`: List of delays of each current pulse in milliseconds

* :const:`LCNhm_configurationfile.CURRENT_AMPLITUDES`: List of amplitudes of each current pulse in nanoampers (nA)

* :const:`LCNhm_configurationfile.CURRENT_SECTION`: List of sections of each current pulse

* :const:`LCNhm_configurationfile.CURRENT_LOCATION`: List of location along the defined :const:`LCNhm_configurationfile.CURRENT_SECTION` of each current pulse

Their length must be the same, the number of different current clamps

Click any of the links for further information
"""
IntrinsicFactors = [CELLPROP_INTRINSIC_IONCHS, CELLPROP_INTRINSIC_EXPERIMENT, CELLPROP_INDIVIDUAL]
"""
List: 
List of all intrinsic properties, packed to be a one of the :class:`LCNhm_class.neuron_class`'s inputs

* :const:`LCNhm_configurationfile.CELLPROP_INTRINSIC_IONCHS`: Ion channels to include in the cell

* :const:`LCNhm_configurationfile.CELLPROP_INTRINSIC_EXPERIMENT`: Additional factor multiplying the following maximum conductances and axial resistance

* :const:`LCNhm_configurationfile.CELLPROP_INDIVIDUAL`: Number of the intrinsic genetic-algorithm `individual`

Click any of the links for further information
"""
SynapticFactors = [CELLPROP_SYNAPTIC_INPUTS, CELLPROP_SYNAPTIC_EXPERIMENT,
				   SIMPROP_THETA_MODE, SIMPROP_THETA_PERIOD, SIMPROP_START_TIME, SIMPROP_END_TIME, SIMPROP_DT]
"""
List:
List of all synaptic properties, packed to be a one of the :class:`LCNhm_class.neuron_class`'s inputs
* :const:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_INPUTS` : Synaptic inputs to include in the cell

* :const:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_EXPERIMENT` : Additional factor multiplying the following maximum conductances

* :const:`LCNhm_configurationfile.SIMPROP_THETA_MODE` : Set theta (``True``) or not (``False``)

* :const:`LCNhm_configurationfile.SIMPROP_THETA_PERIOD` : Theta period in milliseconds

* :const:`LCNhm_configurationfile.SIMPROP_START_TIME` : Lapse of time before starting the simulation in milliseconds

* :const:`LCNhm_configurationfile.SIMPROP_END_TIME` : Total duration of the simulation in milliseconds (:const:`LCNhm_configurationfile.SIMPROP_START_TIME` + :const:`LCNhm_configurationfile.SIMPROP_SIM_TIME`)

Click any of the links for further information
"""
# Make neuron
Pyramidal = neuron_class(MorphoName = CELLPROP_MORPHOLOGY,
						 IntrinsicFactors = IntrinsicFactors,
						 SynapticFactors = SynapticFactors,
						 CurrentFactors = CurrentFactors,
						 DirLocation = DIR_LOCATION ) 
"""
Object: Pyramidal from :class:`LCNhm_class.neuron_class`"""
print '  %3d:%.2d ... pyramidal cell made'%((time.time()-TimeCountStart)/60.,round((time.time()-TimeCountStart)%60))


# 	--------------
# 	Set recordings
# 	--------------
Recordings = recordings(Pyramidal, RECORDING_MAGNITUDE, RECORDING_SECTION, RECORDING_LOCATION)
"""Dictionary: Ouput from :func:`LCNhm_function.recordings`
with recordings of the :class:`LCNhm_class.neuron_class` Pyramidal"""
print '  %3d:%.2d ... recorings set'%((time.time()-TimeCountStart)/60.,round((time.time()-TimeCountStart)%60))


# 	--------------
# 	Run simulation
# 	--------------
# Set simulation parameters
h.tstop = SIMPROP_END_TIME
h.celsius = SIMPROP_TEMPERATURE
# Run
h.run()
print '  %3d:%.2d ... end of simulation'%((time.time()-TimeCountStart)/60.,round((time.time()-TimeCountStart)%60))


# 	------------
# 	Save results
# 	------------
# Detect and save spiking times, recordings and parameters
save_spiking_times(Recordings, FolderName)
save_recordings(Recordings, FolderName)
save_parameters(FolderName)
print '  %3d:%.2d ... results saved'%((time.time()-TimeCountStart)/60.,round((time.time()-TimeCountStart)%60))
