"""
The only script thought to be modified. 
It contains all the simulation parameters:

It's convenient to read the specifications in the documentation, 
because error may arise if the format is not the adequate. 

For instance, if you want to read more about the possible magnitudes 
that can be recorded, you can to go :ref:`search`, and type "RECORDING_MAGNITUDE".
"""

# MAIN DIRECTORY
# --------------
DIR_LOCATION = '/home/andrea/Projects/HippoModel/LCNhippomodel'
OPT_FOLDER_NAME = 'test'

# SIMULATION PROPERTIES
# ---------------------
SIMPROP_THETA_MODE = True
SIMPROP_THETA_PERIOD = 166.
SIMPROP_START_TIME = 0
SIMPROP_SIM_TIME = 10. * SIMPROP_THETA_PERIOD
SIMPROP_END_TIME = SIMPROP_START_TIME + SIMPROP_SIM_TIME
SIMPROP_DT = 0.025
SIMPROP_TEMPERATURE = 34.

# CELL CHARACTERISTICS
# --------------------
CELLPROP_MORPHOLOGY = 'n128'
CELLPROP_INTRINSIC = 2
CELLPROP_SYNAPTIC = 0
CELLPROP_INTRINSIC_IONCHS = ['iNas','iA','iAHPs','iC','iCaL','iCaT','iKDR','iM','iHCN','iL']
CELLPROP_INTRINSIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ]
CELLPROP_SYNAPTIC_INPUTS = ['CA3','CA2','EC3','EC2','Axo','Bis','CCK','Ivy','NGF','OLM','PV','SCA']
CELLPROP_SYNAPTIC_EXPERIMENT = [ 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1. ]

# CURRENT OPTIONS
# ---------------
CURRENT_DURATION = [0.]
CURRENT_DELAY = [0.]
CURRENT_AMPLITUDES = [0.]
CURRENT_SECTION = ['SomaList0']
CURRENT_LOCATION = [0.]

# RECORDING OPTIONS
# -----------------
RECORDING_MAGNITUDE = ['Vmem','Pos']
RECORDING_SECTION = ['SomaList0','ApicList0','ApicList0']
RECORDING_LOCATION = [0.0,0.2,0.6]