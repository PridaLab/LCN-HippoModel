# LCN-HippoModel

One of the most studied brain rhythm is the **theta rhythm**. Theta (4-12 Hz) is considered
the “on-line” status of the hippocampus, a highly rhythmic activity that acts as a global 
synchronizer mechanism for encoding and information processing. During a theta cycle different
neuronal populations of the hippocampus fire sequentially: each one at a distinct preferred phase. 
Moreover, it has been recently discovered that pyramidal neurons of the CA1 hippocampal region, 
generally considered a homogeneous population, can be classified into **deep** and **superficial**:
not only each cell type shows a different preferred phase along theta cycles, but also a different
response in other hippocampal rhythm, sharp wave-ripples, an oscillation associated to the 
consolidation of memory. Understanding how this coordination holds up hippocampal
functions as spatial navigation and memory, is a major question in the field.


The LCN-HippoModel is a **biophysically realistic model** of **CA1 pyramidal cells** aimed to get
novel insights on firing dynamics in deep and superficial populations during the theta rhythm, and
the role of the differential contribution of both the excitatory and inhibitory synaptic inputs, 
and the biophysical intrinsic properties.


The model includes known **excitatory and inhibitory inputs**, using morphologically reconstructions
from NeuroMorpho (http://neuromorpho.org/) (a public database), and a precise distribution of the 
main ion channels via the Hodking-Huxley multi-compartment formalism in the NEURON+Python 
(https://www.neuron.yale.edu/neuron/) platform.

To provide **diversity** among the pyramidal cells, we generated several sets of intrinsic properties
and, called **individuals**, through a genetic algorithm (GA). The GA allowed us to identify a range
of ionic conductances (called genes in the genetic algorithm terminology) that target experimental
values in each given morphology. So our *individuals* are some of those *set of genes* that targeted
experimental constraints. This procedure was repeated with the synaptic properties, so we ended up with
several **intrinsic and synaptic individuals** that would give us heterogenous responses.


The LCN-HippoModel offers an easy way to test hypotheses about the role of different synaptic and intrinsic
factors in the firing dynamics of CA1 pyramidal cells by simply changing experimental parameters defined in the
LCNhm_configurationfile.py.