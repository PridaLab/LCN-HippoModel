"""
The only script thought to be modified. It contains all the simulation parameters:

* **DIR_LOCATION** Directory location where all the model is (eg.: */home/acnavasolive/Projects/LCNhippomodel*)

* **OPT_FOLDER_NAME** Optional name for the folder in which results are going to be saved (eg.: *testing_decreasing_CA3*)

* **SIMPROP_X** Simulation properties, such the ``SIMPROP_START_TIME`` or the ``SIMPROP_TEMPERATURE``.

* **CELLPROP_X** Artificial cell properties, such ``CELLPROP_MORPHOLOGY`` or `CELLPROP_SYNAPTIC_INPUTS``

* **CURRENT_X** Current clamp properties, such the ``CURRENT_DURATION`` or ``CURRENT_AMPLITUDE``.

* **RECORDING_X** Physical magnitudes to be measured and recorded, such ``RECORDING_MAGNITUDE`` or ``RECORDING_LOCATION``.

The script in the in format:

.. code-block:: python
	:linenos:

	DIR_LOCATION = '/home/Projects/LCNhippomodel'
	OPT_FOLDER_NAME = 'test'
	SIMPROP_THETA_MODE = True
	SIMPROP_THETA_PERIOD = 166.
	# ...
	CELLPROP_INTRINSIC = 0
	CELLPROP_INTRINSIC_IONCHS = ['iNas','iA','iAHPs','iC','iCaL','iCaT','iKDR','iM','iHCN','iL']
	# ...

so in order to change the parameters it's just writing over the default values.

It's convenient to read the specifications in the documentation, because error may arise if the format is not the adequate. For instance, if you want to read more about the possible magnitudes that can be recorded, you can to go :ref:`search`, and type "RECORDING_MAGNITUDE".


Parameters and descriptions
---------------------------

"""

# MAIN DIRECTORY
# --------------
DIR_LOCATION = '/home/andrea/Projects/HippoModel/LCNhippomodel'
"""
String: Precise location of the main directory.

	Eg.: /home/acnavasolive/Projects/LCNhippomodel
"""
OPT_FOLDER_NAME = 'test'
"""
OPT_FOLDER_NAME: String
    Optional folder name: additional description of the simulation.

    Defined in :const:`LCNhm_configurationfile.OPT_FOLDER_NAME`
    
    Eg.: *first_test*, *high_CA3*, *inhibition_suppression*
"""

# SIMULATION PROPERTIES
# ---------------------
SIMPROP_THETA_MODE = True
""" Boolean: Statement for setting theta state"""
SIMPROP_THETA_PERIOD = 166.
""" Float: Theta period in milliseconds"""
SIMPROP_START_TIME = 0
""" Float: Lapse of time before starting the simulation in milliseconds. 

.. warning::
	SIMPROP_SIM_TIME must be greater or equal to 1.5 * SIMPROP_THETA_PERIOD,
	if SIMPROP_THETA_MODE is  True

"""
SIMPROP_SIM_TIME = 10. * SIMPROP_THETA_PERIOD
""" Float: Duration of the simulation in milliseconds"""
SIMPROP_END_TIME = SIMPROP_START_TIME + SIMPROP_SIM_TIME
""" Float: Total duration of the simulation in milliseconds (``SIMPROP_START_TIME`` + ``SIMPROP_SIM_TIME``)"""
SIMPROP_DT = 0.025
""" Float: Time step in milliseconds

Default is 0.025 ms"""
SIMPROP_TEMPERATURE = 34.
""" Float: Temperature of simulation (celsius)"""


# CELL CHARACTERISTICS
# --------------------
CELLPROP_MORPHOLOGY = 'n128'
"""
String: Morphology name, choose between `n128`, `sup1` or `n409`
"""
CELLPROP_INTRINSIC = 2
"""
Int: Number of the intrinsic-genetic algorithm `individual`
	
	Choose between 0, 1, 2, ..., 12

Individuals are defined as factors multiplying the following maximum conductances and axial resistance:

	gNas,  gA,  gAHPs,  gC,  gKDR,  gM,  gCaL,  gCaT,  gHCN,  gL (all in S/cm2),  and Ra (Ohm*cm)

.. note::
	Not all combinations of intrinsic and synaptic `individuals` are validated
	Please use one the following combinations:

	* `CELLPROP_INTRINSIC = 1`: `CELLPROP_SYNAPTIC` from {1}

	* `CELLPROP_INTRINSIC = 2`: `CELLPROP_SYNAPTIC` from {1, 3}

	* `CELLPROP_INTRINSIC = 3`: `CELLPROP_SYNAPTIC` from {0, 1}

	* `CELLPROP_INTRINSIC = 4`: `CELLPROP_SYNAPTIC` from {1, 3}

	* `CELLPROP_INTRINSIC = 5`: `CELLPROP_SYNAPTIC` from {3}

	* `CELLPROP_INTRINSIC = 6`: `CELLPROP_SYNAPTIC` from {0}

	* `CELLPROP_INTRINSIC = 7`: `CELLPROP_SYNAPTIC` from {0, 3}

	* `CELLPROP_INTRINSIC = 8`: `CELLPROP_SYNAPTIC` from {3}

	* `CELLPROP_INTRINSIC = 9`: `CELLPROP_SYNAPTIC` from {0, 1}

	* `CELLPROP_INTRINSIC = 10`: `CELLPROP_SYNAPTIC` from {0, 1, 3}

	* `CELLPROP_INTRINSIC = 11`: `CELLPROP_SYNAPTIC` from {1}

	* `CELLPROP_INTRINSIC = 12`: `CELLPROP_SYNAPTIC` from {0, 1}

"""
CELLPROP_SYNAPTIC = 0
"""
String: Number of the synaptic-genetic algorithm `individual`
	
	Choose between 0, 1, 2, 3

.. note::
	Not all combinations of synaptic and intrinsic `individuals` are validated
	Please use one the following combinations:

	* `CELLPROP_SYNAPTIC = 0`: `CELLPROP_INTRINSIC` from {3, 6, 7, 9, 10, 12}

	* `CELLPROP_SYNAPTIC = 1`: `CELLPROP_INTRINSIC` from {1, 2, 3, 4, 9, 10, 11, 12}

	* `CELLPROP_SYNAPTIC = 3`: `CELLPROP_INTRINSIC` from {2, 4, 5, 7, 8, 10}
"""
CELLPROP_INTRINSIC_IONCHS = ['iNas','iA','iAHPs','iC','iCaL','iCaT','iKDR','iM','iHCN','iL']
"""List of Strings: Ion channels to include in the cell.

Default (All) - iA, iAHPs, iC, iCaL, iCaT, iHCN, iKDR, iL, iM, iNas

* iA - A-type K+ clannels (Cutsuridis2015)

* iAHPs - CA2+ dependent slow AHP K+ conductance (Cutsuridis2015)

* iC - Short-duration [Ca]- and voltage-dependent current (Cutsuridis2015)
	
* iCaL - L-type calcium channel with low threshold for activation (Muellner2015)
	
* iCaT - Somatic L-type calcium channel with low threshold for activation (Muellner2015)
	
* iHCN - Hyperpolarization-activated cyclic nucleotide-gated channels (Sinha2015)
	
* iKDR - Delay rectifier current (Cutsuridis2015)
	
* iL - Leak channels (default)
	
* iM - Slowly activating voltage-dependent potassium current (Migliore2006)
	
* iNas - Slow sodium channels (Jaslove1992)
"""
CELLPROP_INTRINSIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ]
"""
List of Floats: Additional factor multiplying the following maximum conductances and axial resistance
	
	gNas,  gA,  gAHPs,  gC,  gKDR,  gM,  gCaL,  gCaT,  gHCN,  gL (all in S/cm2),  and Ra (Ohm*cm)

Default, all factors equal 1. To experiment increasing/decreasing intrinsic properties, change the appropiate element of the list

	Eg.: Decrease Nas conductance to zero

	>>> CELLPROP_INTRINSIC_EXPERIMENT = [ 0., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ]

	Eg.: Increase axial resistance and HCN conductance x2

	>>> CELLPROP_INTRINSIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 2., 1., 2. ]
"""
CELLPROP_SYNAPTIC_INPUTS = ['CA3','CA2','EC3','EC2','Axo','Bis','CCK','Ivy','NGF','OLM','PV','SCA']
"""
List of Strings: Synaptic inputs to include in the cell (Bezaire2016)

Default (All) - CA3,  CA2,  EC3,  EC2,  Axo,  Bis,  CCK,  Ivy,  NGF,  OLM,  PV,  SCA
"""
CELLPROP_SYNAPTIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ]
"""
List of Floats: Additional factor multiplying the following maximum conductances
	
	gCA3,  gCA2,  gEC3,  gEC2,  gAxo,  gBis,  gCCK,  gIvy,  gNGF,  gOLM,  gPV,  gSCA (uS)

Default, all factors equal 1. To experiment increasing/decreasing synaptic properties, change the appropiate element of the list.

	Eg.: Superficial pyramidal

	>>> CELLPROP_SYNAPTIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 0.3, 1. ]
	
	Eg.: Deep pyramidal

	>>> CELLPROP_SYNAPTIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 0.3, 1., 1., 1., 1. ]

	Eg.: Just excitation

	>>> CELLPROP_SYNAPTIC_EXPERIMENT = [ 1., 1., 1., 1., 0., 0., 0., 0., 0., 0., 0. ]
"""

# CURRENT OPTIONS
# ---------------
CURRENT_DURATION = [0.]
"""List of Floats: Durations of each current pulse in milliseconds

	Eg.: Three pulses during 50, 100 and 150ms

	>>> CURRENT_DURATION = [50., 100., 150.]

	Eg.: No pulse

	>>> CURRENT_DURATION = [0.]

Length of :const:`CURRENT_DURATION` must be the same to :const:`CURRENT_AMPLITUDES`, :const:`CURRENT_DELAY`, :const:`CURRENT_SECTION` and :const:`CURRENT_LOCATION`
"""
CURRENT_DELAY = [0.]
"""
List of Floats: Delays of each current pulse in milliseconds

	Eg.: Three pulses starting after 500ms of the simulation start, separated 100ms

	>>> CURRENT_DELAY = [500., 600., 700.]

	Eg.: No pulse

	>>> CURRENT_DELAY = [0.]

Length of :const:`CURRENT_DELAY` must be the same to :const:`CURRENT_AMPLITUDES`, :const:`CURRENT_DURATION`, :const:`CURRENT_SECTION` and :const:`CURRENT_LOCATION`
"""
CURRENT_AMPLITUDES = [0.]
"""
List of Floats: Amplitudes of each current pulse in nanoampers (nA)

	Eg.: Three pulses of 0.1, 0.2 and 0.3 nA

	>>> CURRENT_AMPLITUDES = [0.1, 0.2, 0.3]

	Eg.: No pulse

	>>> CURRENT_AMPLITUDES = [0.]

Length of :const:`CURRENT_AMPLITUDES` must be the same to :const:`CURRENT_DELAY`, :const:`CURRENT_DURATION`, :const:`CURRENT_SECTION` and :const:`CURRENT_LOCATION`
"""
CURRENT_SECTION = ['SomaList0']
"""
List of ``str``: Sections of each current pulse. Choose between:

* n128
	* Soma: SomaList0, ..., SomaList3
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList30
	* Apical dendrites: ApicList0, ..., ApicList67
	* Somato-apical dendrites: ApicList0

* sup1
	* Soma: SomaList0
	* Axon: AxonList0, AxonList1
	* Basal dendrites: DendList0, ..., DendList46
	* Apical dendrites: ApicList0, ..., ApicList48
	* Somato-apical dendrites: ApicList4

* n409
	* Soma: SomaList0, SomaList1
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList41
	* Apical dendrites: ApicList0, ..., ApicList40
	* Somato-apical dendrites: ApicList0

* n127
	* Soma: SomaList0, ..., SomaList4
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList22
	* Apical dendrites: ApicList0, ..., ApicList44
	* Somato-apical dendrites: ApicList0



	Eg.: Three pulses of 0.1nA, in the middle of ``SomaList0``, ``ApicList20`` and ``DendList10``, simultaneouly

	>>> CURRENT_SECTION = [SomaList0, ApicList20, DendList10]
	>>> CURRENT_LOCATION = [0.5, 0.5, 0.5]
	>>> CURRENT_AMPLITUDES = [0.1, 0.1, 0.1]
	>>> CURRENT_DELAY = [0., 0., 0.]
	>>> CURRENT_DURATION = [100, 100, 100]

	Eg.: No pulse (anything)

	>>> CURRENT_SECTION = [SomaList0]

Length of :const:`CURRENT_SECTION` must be the same to :const:`CURRENT_DELAY`, :const:`CURRENT_DURATION`, :const:`CURRENT_AMPLITUDES` and :const:`CURRENT_LOCATION`
"""
CURRENT_LOCATION = [0.]
"""
List of ``str``: Location along the defined :const:`CURRENT_SECTION` of each current pulse. 

Location is defined between 0 and 1, being 0 the beginning on the ``Section``, and 1 the end.

	Eg.: Two pulses of 0.1nA at the beginning and the middle of ``ApicList0``

	>>> CURRENT_LOCATION = [0.0, 0.5]
	>>> CURRENT_SECTION = [ApicList0, ApicList0]
	>>> CURRENT_AMPLITUDES = [0.1, 0.1]
	>>> CURRENT_DELAY = [0., 0.]
	>>> CURRENT_DURATION = [100, 100]

	Eg.: No pulse (anything)

	>>> CURRENT_LOCATION = [0.]

Length of :const:`CURRENT_LOCATION` must be the same to :const:`CURRENT_DELAY`, :const:`CURRENT_DURATION`, :const:`CURRENT_SECTION` and :const:`CURRENT_AMPLITUDES`
"""

# RECORDING OPTIONS
# -----------------
RECORDING_MAGNITUDE = ['Vmem','Pos']
"""
List of ``str``: Physical magnitudes to be measured and recorded.
Options are:

* Time: Simulation time (ms)

* Vmem: Membrane potential (mV)

* Imem: Membrane current (mA/cm2)

* Pos: Position by three spatial coordinates (us), being *z-axis* the radial axis

	Eg.: Save only membrane potential at soma and some points at the somato-apical trunk

	>>> RECORDING_MAGNITUDE = ['Vmem']
	>>> RECORDING_SITE = ['SomaList0','SomApicList1','SomApicList3','SomApicList5']

	Eg.: Save membrane potential, current and position at soma

	>>> RECORDING_MAGNITUDE = ['Vmem','Imem','Pos']
	>>> RECORDING_SITE = ['SomaList0']
"""
RECORDING_SECTION = ['SomaList0','ApicList0','ApicList0']
"""
List of ``str``: Sites/Sections to be measured and recorded.
Options are:

* n128
	* Soma: SomaList0, ..., SomaList3
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList30
	* Apical dendrites: ApicList0, ..., ApicList67
	* Somato-apical dendrites: ApicList0

* sup1
	* Soma: SomaList0
	* Axon: AxonList0, AxonList1
	* Basal dendrites: DendList0, ..., DendList46
	* Apical dendrites: ApicList0, ..., ApicList48
	* Somato-apical dendrites: ApicList4

* n409
	* Soma: SomaList0, SomaList1
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList41
	* Apical dendrites: ApicList0, ..., ApicList40
	* Somato-apical dendrites: ApicList0

* n127
	* Soma: SomaList0, ..., SomaList4
	* Axon: AxonList0
	* Basal dendrites: DendList0, ..., DendList22
	* Apical dendrites: ApicList0, ..., ApicList44
	* Somato-apical dendrites: ApicList0

	Eg.: Save only membrane potential at soma and somato-apical trunk

	>>> RECORDING_MAGNITUDE = ['Vmem']
	>>> RECORDING_SITE = ['SomaList0','ApicList0']

	Eg.: Save membrane potential, current and position at soma

	>>> RECORDING_MAGNITUDE = ['Vmem','Imem','Pos']
	>>> RECORDING_SECTION = ['SomaList0']

Lengths of :const:`RECORDING_SECTION` and :const:`RECORDING_LOCATION` must be the same
"""
RECORDING_LOCATION = [0.0,0.2,0.6]
"""
List of Floats: Location along the defined :const:`RECORDING_SECTION` to be measured and recorded

Location is defined between 0 and 1, being 0 the beginning on the ``Section``, and 1 the end.

	Eg.: Save only membrane potential at soma and at the beginning, middle and end of somato-apical trunk

	>>> RECORDING_MAGNITUDE = ['Vmem']
	>>> RECORDING_SECTION = ['SomaList0','ApicList0','ApicList0','ApicList0']
	>>> RECORDING_LOCATION = [     0.0,        0.0,        0.5,        1.0  ]

	Eg.: Save membrane potential, current and position at some basal dendrites

	>>> RECORDING_MAGNITUDE = ['Vmem','Imem','Pos']
	>>> RECORDING_SECTION = ['DendList0','DendList3','DendList5','DendList10']
	>>> RECORDING_LOCATION = [     0.3,        0.3,        0.3,         0.3  ]

Lengths of :const:`RECORDING_SECTION` and :const:`RECORDING_LOCATION` must be the same
"""
