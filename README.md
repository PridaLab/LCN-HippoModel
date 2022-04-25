# LCN-HippoModel

## Place field branch

As an addition to a constant theta state, we have added this branch in which a modulation is added to simulate crossing a place field. For that, we added a place field (PF) asymmetric modulation to CA3, axo-axonic, PV-bc and CCK-bc synaptic conductances. The center of the place field is assumed to be at 72% of place field width (Harvey, 2009). In order to inform the extent into which CA3, Axo, PVbc and CCKbc are modulated, a new input to the neuron class `neuron_class` is added, called `DGiDGe`, which is a two-length array `[ ΔG(inhibition), ΔG(excitation) ]` with:

	ΔG = ((max G inside PF) - (constant G outside PF)) / (constant G outside PF)

so if: 
	* `ΔG = 0`: No modulation
	* `ΔG = 1`: Conductance doubles inside place field
	* `ΔG = -1`: Conductance decreases to 0 inside place field

The current way to use the class would be:

		Pyramidal = neuron_class_PFmodulation(MorphoName = CELLPROP_MORPHOLOGY,
								 IntrinsicFactors = IntrinsicFactors,
								 SynapticFactors = SynapticFactors,
								 CurrentFactors = CurrentFactors,
								 DirLocation = DIR_LOCATION,
								 DGiDGe = [DGi, DGe]) 

# Documentation

Full documentation of the code can be found here:

https://acnavasolive.github.io/LCN-HippoModel/


# Full article

You can read [here](https://www.nature.com/articles/s41467-020-15840-6) the article where this model was used.

**Navas-Olive, A.**, Valero, M., Jurado-Parras, T. et al. Multimodal determinants of phase-locked dynamics across deep-superficial hippocampal sublayers during theta oscillations. _Nat Commun_ 11, 2217 (2020). 
