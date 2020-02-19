"""
This script shapes the class from which our neurons will be built. 

Neurons will have the following attributes:

* **MorphoName**: Morphology name, it must be one of the following: `n128`, `sup1`, `n409` or `n127`
* **MorphoData**: Dataset containing all the morphological data
* **SomaList**: List of NEURON somatic ``Sections``
* **AxonList**: List of NEURON axonic ``Sections``
* **DendList**: List of NEURON Basal dendritic ``Sections``
* **ApicList**: List of NEURON apical dendritic ``Sections``
* **SomApicList**: List of somatoapical NEURON apical dendritic ``Sections``
* **SectionsList**: List of all NEURON ``Sections``, the union of :attr:`SomaList`, :attr:`AxonList`, :attr:`DendList`, :attr:`ApicList`
* **TopolDict**: Dictionary with all topological information in the form
* **CurrentObject**: List of NEURON objects containing all current clamps information


Access to all these properties is possible by:

.. code-block:: python

    Eg.:
    
    >>> Pyramidal = neuron_class(MorphoName = CELLPROP_MORPHOLOGY,
                 IntrinsicFactors = IntrinsicFactors,
                 SynapticFactors = SynapticFactors,
                 CurrentFactors = CurrentFactors,
                 DirLocation = DIR_LOCATION ) 

    >>> print Pyramidal.MorphoName

    sup1


The way the neuron is built is following these steps:

1. **Import and set morphological data** from ./LCNhm-neurondata/``CELLPROP_MORPHOLOGY``.swc
    These .swcs cointain information about the 3D structure and width of the cell. 

    You can obtain more morphologies on the public databas `NeuroMorpho <http://neuromorpho.org/>`_

    Classification into four main morphological categories (soma, axon, basal and apical dendrites).

    Every category will be compartimentalized in multiple sections, each of which will have different and specific
    ion channels, synaptic inputs, diameter, spines, etc... For instance, one apical dendrite that starts on the somatoapical trunk,
    will be compartimentalized in 10 sections, where the nearest section to the trunk will be thicker than the furthest.

#. **Set biophysics** (ion channels, membrane capacitance, membrane resistivity, etc...)
    What ion channels to use can be chosen, as well as a factor to increase/decrease the conductance density.

    See :data:`LCNhm_configurationfile.CELLPROP_INTRINSIC_IONCHS` and :data:`LCNhm_configurationfile.CELLPROP_INTRINSIC_EXPERIMENT`: 
    for further details.

#. **Set synapses** (excitatory/inhibitory inputs, place of synaptic boutons, maximum conductance, etc...)
    Synaptic inputs can be chosen, as well as a factor to increase/decrease the maximum conductance.

    See :data:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_INPUTS` and :data:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_EXPERIMENT`: 
    for further details.

#. **Set current clamp** (place, duration, intensity, etc...)
    A square-pulse clamp is set according to parameters written in :ref:`LCNhm-configuration-file`.

    See :data:`LCNhm_configurationfile.CURRENT_DURATION`, :data:`LCNhm_configurationfile.CURRENT_DELAY`, :data:`LCNhm_configurationfile.CURRENT_AMPLITUDES`, :data:`LCNhm_configurationfile.CURRENT_SECTION`, :data:`LCNhm_configurationfile.CURRENT_LOCATION`, for further details.
"""
import sys
import math
import warnings
import numpy as np
import pandas as pd
from neuron import h, nrn, gui
from LCNhm_functions import *


class neuron_class(object):
    """
    Neurons will be built as `neuron_class` objects.

    Parameters
    ----------
        MorphoName: String
            Morphology name, it must be one of the following: `n128`, `sup1`, `n409` or `n127`

            Defined in :const:`LCNhm_configurationfile.CELLPROP_MORPHOLOGY`

        IntrinsicFactors: List
            List of all current properties

            Defined in :const:`LCNhm_main.IntrinsicFactors`

        SynapticFactors: List
            List of all intrinsic properties

            Defined in :const:`LCNhm_main.SynapticFactors`

        CurrentFactors: List
            List of all synaptic properties

            Defined in :const:`LCNhm_main.CurrentFactors`

        DirLocation: String
            Precise location of the main directory

    Returns
    -------
        ``neuron_class`` object: Class object

    Attributes
    ----------
        MorphoName: String
            Morphology name, it must be one of the following: `n128`, `sup1`, `n409` or `n127`

            Defined in :const:`LCNhm_configurationfile.CELLPROP_MORPHOLOGY`

            It is one of the :class:`neuron_class` inputs

        MorphoData: csv
            Dataset containing all the morphological data

            It describes each morphological point with six parameters:
                * ``Type``: soma (1), axon (2), Basal (3) or apical (4) dendrite
                * ``x`` , ``y`` , ``z``: Spatial coordinates (in micrometers)
                * ``d``: diameter (in micrometers)
                * ``IDFather``: Line number of the morphological point to which is connected

            See :func:`make_geometry_dictionary` to see how it is used

            Source file in *../LCNhm-neurondata/ <MorphoName> .swc*

        SomaList: List
            List of NEURON somatic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        AxonList: List
            List of NEURON axonic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        DendList: List
            List of NEURON Basal dendritic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        ApicList: List
            List of NEURON apical dendritic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        SomApicList: List
            List of somatoapical NEURON apical dendritic ``Sections``

            List can be found in *./LCNhm-neurondata/somatoapical_sections.txt*

        SectionsList: List
            List of all NEURON ``Sections``, the union of :attr:`SomaList`, :attr:`AxonList`, :attr:`DendList`, :attr:`ApicList`

        TopolDict: Dictionary
            Dictionary with all topological information in the form

                [ [ SecListName[SecNumber] ,  Who is connected to (IDFather) ,  x ,  y ,  z ] ] , [...], ... ]

            Eg.:

            .. code-block:: python

                #            NEURON object    name           ID of father     x  y  z
                TopolDict = { SomaList[0]: ['SomaList[0]',       -1,          0, 0, 0],
                              SomaList[1]: ['SomaList[1]', ID of SomaList[0], 0, 0, -2],
                              ...
                              ApicList[9]: ['ApicList[9]', ID of ApicList[2], 0, 0, -1000]}

        CurrentObject: List
            List of NEURON objects containing all current clamps information
            
                [[ Clamp object 1, Duration 1, Delay 1, Amplitude 1 ],
                [  Clamp object 2, Duration 2, Delay 2, Amplitude 2 ], [...], ...... ]

            Eg.:

            .. code-block:: python

                #                   Object                 dur   del  amp 
                CurrentObject = [ [ NEURON object at soma, 100,    0,   0],
                                  [ NEURON object at soma, 200,  300,   0],
                                  ...
                                  [ NEURON object at soma, 200, 1200,   0]]

    """

    # Atributes and callings to methods
    def __init__(self, MorphoName, IntrinsicFactors, SynapticFactors, CurrentFactors, DirLocation):

        global DirLoc; DirLoc = DirLocation

        # IMPORTING DATA
        # ==============
        # Name of the morphology
        self.MorphoName = MorphoName
        # Dataset containing all the morphological data
        self.MorphoData = pd.read_csv('%s/LCNhm-neurondata/%s.swc'%(DirLoc,MorphoName),
                                    names = ['Type','x','y','z','d','IDFather'],
                                    delimiter = ' ', index_col = 0)

        # SECTIONS
        # ========
        # Initialization and group in lists of the four main class of sections: soma, axon, Basal (dend) and apical (apic) dendrites
        self.SomaList, self.AxonList, self.DendList, self.ApicList, self.SomApicList = self.init_sections(self.MorphoData)
        # All sections in one list
        self.SectionsList = self.SomaList+self.AxonList+self.DendList+self.ApicList

        # GEOMETRY AND TOPOLOGY
        # =====================
        # TopolDict is the dictionary with the geometry
        self.CellPosition = [0, 0, 0, 0]
        self.TopolDict = self.make_geometry_dictionary(self.CellPosition, self.MorphoData)
        # Making the connections
        self.set_geometry(self.TopolDict)

        # BIOPHYSICS 
        # ==========
        # We define the biophysics
        self.set_intrinsic_properties(IntrinsicFactors)
        
        # SYNAPSES
        # ========
        self.SynapticFactors = SynapticFactors
        #self.SynDict, self.synObj, self.NetConObjects = self.SetSynapses(self.SynapticFactors)
        if len(SynapticFactors[0])>0: self.SynDict, self.SynObjects, self.NetConObjects = self.set_synaptic_properties(SynapticFactors)

        # ICLAMP
        # ======
        self.CurrentObject = self.set_current_clamp(CurrentFactors) 


    def init_sections(self, MorphoData):
        """
        Initialization and group in lists of the four main class of sections: soma, axon, Basal (dend) and apical (apic) dendrites

        Parameters
        ----------
        MorphoData: csv
            :class:`neuron_class`'s :attr:`MorphoData` attribute: Dataset containing all the morphological data

            It describes each morphological point with six parameters:
                * ``Type``: soma (1), axon (2), Basal (3) or apical (4) dendrite
                * ``x`` , ``y`` , ``z``: Spatial coordinates (in micrometers)
                * ``d``: diameter (in micrometers)
                * ``IDFather``: Line number of the morphological point to which is connected

            See :func:`make_geometry_dictionary` to see how it is used

            Source file in *../LCNhm-neurondata/ <MorphoName> .swc*

        Returns
        -------
        SomaList: List
            class:`neuron_class`'s :attr:`SomaList` attribute: List of NEURON somatic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        AxonList: List
            class:`neuron_class`'s :attr:`AxonList` attribute: List of NEURON axonic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        DendList: List
            class:`neuron_class`'s :attr:`DendList` attribute: List of NEURON Basal dendritic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        ApicList: List
            class:`neuron_class`'s :attr:`ApicList` attribute: List of NEURON apical dendritic ``Sections``

            The number of ``Sections`` is derived from :attr:`MorphoData`

        SomApicList: List
            class:`neuron_class`'s :attr:`SomApicList` attribute: List of somatoapical NEURON apical dendritic ``Sections``.

            List can be found in `./LCNhm-neurondata/somatoapical_sections.txt`.
        """
        # We count how many lines have Type 1 (SomaList)
        numsoma = sum([1 for i in MorphoData.index if (MorphoData['Type'][i] == 1) and ((MorphoData['IDFather'][i] != i-1) or (MorphoData['Type'][MorphoData['IDFather'][i]] != MorphoData['Type'][i]))])
        # We count how many lines have Type 2 (AxonList)
        numaxon = sum([1 for i in MorphoData.index if (MorphoData['Type'][i] == 2) and ((MorphoData['IDFather'][i] != i-1) or (MorphoData['Type'][MorphoData['IDFather'][i]] != MorphoData['Type'][i]))])
        # We count how many lines have Type 3 (DendList)
        numdend = sum([1 for i in MorphoData.index if (MorphoData['Type'][i] == 3) and ((MorphoData['IDFather'][i] != i-1) or (MorphoData['Type'][MorphoData['IDFather'][i]] != MorphoData['Type'][i]))])
        # We count how many lines have Type 4 (ApicList)
        numapic = sum([1 for i in MorphoData.index if (MorphoData['Type'][i] == 4) and ((MorphoData['IDFather'][i] != i-1) or (MorphoData['Type'][MorphoData['IDFather'][i]] != MorphoData['Type'][i]))])

        # We create four section lists: soma, axon, Basal and apical dendrites
        SomaList = [h.Section(name = 'soma', cell = self) for i in range(numsoma)]
        AxonList = [h.Section(name = 'axon', cell = self) for i in range(numaxon)]
        DendList = [h.Section(name = 'dend', cell = self) for i in range(numdend)]
        ApicList = [h.Section(name = 'apic', cell = self) for i in range(numapic)]

        # We create a complmentary section list: somato-apical dendrite
        with open(DirLoc+'/LCNhm-neurondata/somatoapical_sections.txt') as file:
            SomApicList = [ApicList[[int(line.split()[1]) for line in file if (line.split()[0] == self.MorphoName)][0]]]

        # Return
        return SomaList, AxonList, DendList, ApicList, SomApicList

    def make_geometry_dictionary(self, CellPosition, MorphoData):
        """
        Make a legible dictionary from the swc ``MorphoData`` database

        Parameters
        ----------
        CellPosition: List of Floats
            List of four elements

            * x-, y-, z-: spatial position of the soma (in micrometers)

            * angle: self-rotation around its main axis (in degrees)

        MorphoData: csv
            :class:`neuron_class`'s :attr:`MorphoData` attribute: Dataset containing all the morphological data

            It describes each morphological point with six parameters:
                * ``Type``: soma (1), axon (2), Basal (3) or apical (4) dendrite
                * ``x`` , ``y`` , ``z``: Spatial coordinates (in micrometers)
                * ``d``: diameter (in micrometers)
                * ``IDFather``: Line number of the morphological point to which is connected

            See :func:`make_geometry_dictionary` to see how it is used

            Source file in *../LCNhm-neurondata/ <MorphoName> .swc*

        Returns
        -------
        TopolDict: Dictionary
            :class:`neuron_class`'s :attr:`TopolDict` attribute: Dictionary with all topological information in the form

                [ [ SecListName[SecNumber] ,  Who is connected to (IDFather) ,  x ,  y ,  z ] ] , [...], ... ]

            Eg.:

            .. code-block:: python

                #            NEURON object    name           ID of father     x  y  z
                TopolDict = { SomaList[0]: ['SomaList[0]',       -1,          0, 0, 0]
                              SomaList[1]: ['SomaList[1]', ID of SomaList[0], 0, 0, -2],
                              ...
                              ApicList[9]: ['ApicList[9]', ID of ApicList[2], 0, 0, -1000]}

        """


        # Unpack x,y,z,rotation parameters from CellPosition
        X, Y, Z, SelfRotation = CellPosition

        class Point(object):
            """
            Each point (compartment) in the neuron will be a ``Point`` class

            Attributes
            ----------
            ID: Int
                Line number in the MorphoData csv

            Type: Int
                NEURON ``Section`` Type: 1 (soma), 12 (axon), 3 (dend), 4 (apic)

            coords:
                Coordinates from the MorphoData csv, transformed by CellPosition

            IDFather:
                Line number in the MorphoData csv from which this compartment starts
            """
            def __init__(self, ID, Type, coords, sect_num, IDFather, isFather):
                self.ID = ID
                # SomaList, AxonList, Basal Dendrite, Apical Dendrite
                self.Type = Type
                # Rotation
                self.coords = np.array([math.cos(SelfRotation)*coords[0]-math.sin(SelfRotation)*coords[1],
                                        math.sin(SelfRotation)*coords[0]+math.cos(SelfRotation)*coords[1],
                                        coords[2],
                                        coords[3]])
                # Translation
                self.coords += np.array([X, Y, Z, 0])
                self.sect_num = sect_num
                self.IDFather = IDFather
                self.isFather = isFather


        points = []
        numsoma = -1
        numaxon = -1
        numdend = -1
        numapic = -1

        # List of all points
        for ID in MorphoData.index:
            # Type and [x,y,z,diam]
            Type = MorphoData['Type'][ID]
            coords = [MorphoData['x'][ID] , MorphoData['y'][ID] , MorphoData['z'][ID] , MorphoData['d'][ID]]
            IDFather = MorphoData['IDFather'][ID]
            # Is it a father?
            isFather = True if (MorphoData['IDFather'][ID] != (ID-1)) or (MorphoData['Type'][MorphoData['IDFather'][ID]] != Type) else False
            # If it is a father, we include in a new subsection of the section
            if Type == 1: numsoma += 1*isFather; sect_num = numsoma # SomaList
            if Type == 2: numaxon += 1*isFather; sect_num = numaxon # AxonList
            if Type == 3: numdend += 1*isFather; sect_num = numdend # Basal Dendrites
            if Type == 4: numapic += 1*isFather; sect_num = numapic # Apical Dendrites
            # Add point to list
            points.append(Point(ID, Type, coords, sect_num, IDFather, isFather ))

        # Now we make the dictionary:
        TopolDict = {}
        # For each point:
        for point in points:
            # What section Type is it?
            if   point.Type == 1: sect = self.SomaList[point.sect_num]
            elif point.Type == 2: sect = self.AxonList[point.sect_num]
            elif point.Type == 3: sect = self.DendList[point.sect_num]
            elif point.Type == 4: sect = self.ApicList[point.sect_num]
            # If it is a father, we create a new key in the dictionary
            if point.isFather:
                # Point's section name
                if   point.Type == 1: TopolDict[sect] = ['SomaList[%d]'%point.sect_num]
                elif point.Type == 2: TopolDict[sect] = ['AxonList[%d]'%point.sect_num]
                elif point.Type == 3: TopolDict[sect] = ['DendList[%d]'%point.sect_num]
                elif point.Type == 4: TopolDict[sect] = ['ApicList[%d]'%point.sect_num]
                # The first line of the SWC has a "-1" IDFather
                if point.IDFather == -1:
                    grandfather     = False
                    type_grandfather = False
                else:
                    # Where does it come from? As it is a father, grandfather is its father
                    ID_grandfather      = point.IDFather-1 # "-1" because starts counting on 1
                    type_grandfather    = points[ID_grandfather].Type
                    num_grandfather     = points[ID_grandfather].sect_num
                # Grandfather section
                if   type_grandfather == False: TopolDict[sect].append(self.SomaList[point.sect_num])
                elif type_grandfather == 1: 
                    TopolDict[sect].append(self.SomaList[num_grandfather])
                    TopolDict[sect].append([pt.ID for pt in points if (point.ID == pt.IDFather)][0])
                elif type_grandfather == 2: 
                    TopolDict[sect].append(self.AxonList[num_grandfather])
                    TopolDict[sect].append([pt.ID for pt in points if (point.ID == pt.IDFather)][0])
                elif type_grandfather == 3: 
                    TopolDict[sect].append(self.DendList[num_grandfather])
                    TopolDict[sect].append([pt.ID for pt in points if (point.ID == pt.IDFather)][0])
                elif type_grandfather == 4: 
                    TopolDict[sect].append(self.ApicList[num_grandfather])
                    TopolDict[sect].append([pt.ID for pt in points if (point.ID == pt.IDFather)][0])
                # Coordinates of the Point
                TopolDict[sect].append(point.coords)
            else:
                TopolDict[sect].append(point.coords)

        # Once the dictionary is done, we erase the inapropiate dendrites listed in ./LCNhm-neurondata/dendrites_to_erase.txt
        with open(DirLoc+'/LCNhm-neurondata/dendrites_to_erase.txt' ) as file:
            Dend2Erase = [map(int, line.split()[1:]) for line in file if line.split()[0] == self.MorphoName][0]
        self.DendList = [self.DendList[ii] for ii in range(len(self.DendList)) if ii not in Dend2Erase]

        return TopolDict

    def set_geometry(self, TopolDict):
        """ 
        Set the topology and geometry from the :attr:`TopolDict` dictionary.

        * First, the beginning of each ``Section`` is connected to the end of its "father" ``Section``, defined in the :attr:`MorphoData` csv

        * Second, it sets each compartment position and diameter

        Parameters
        ----------
        TopolDict: Dictionary
            :class:`neuron_class`'s :attr:`TopolDict` attribute: Dictionary with all topological information in the form

                [ [ SecListName[SecNumber] ,  Who is connected to (IDFather) ,  x ,  y ,  z ] ] , [...], ... ]

            Eg.:

            .. code-block:: python

                #            NEURON object    name           ID of father     x  y  z
                TopolDict = { SomaList[0]: ['SomaList[0]',       -1,          0, 0, 0]
                              SomaList[1]: ['SomaList[1]', ID of SomaList[0], 0, 0, -2],
                              ...
                              ApicList[9]: ['ApicList[9]', ID of ApicList[2], 0, 0, -1000]}
        """
        for key in TopolDict.keys():
            # The first section (SomaList[0]) doesn't need to be connected
            if key != self.SomaList[0]:
                key.connect(TopolDict[key][1])

        # Once all sections are connected, we define its position and diameter:
        h.pop_section()
        for key in TopolDict.keys():
            if key in self.SomApicList:
                numsegs = 15
                key.push()
                h.pt3dclear()
                pointList = []
                for i in range(3, len(TopolDict[key])):
                    # If it's a list, then add that 3d point
                    pointList.append([TopolDict[key][i][0], TopolDict[key][i][1], TopolDict[key][i][2], TopolDict[key][i][3]])
                ID_compartment = np.array([float(i)/numsegs for i in range(numsegs) ])*len(pointList) 
                for idc in ID_compartment:
                    xi, yi, zi, di = pointList[int(idc)]
                    h.pt3dadd(xi, yi, zi, di)
                h.pop_section()
                key.nseg = numsegs
                continue
            else:
                key.push()
                h.pt3dclear()
                pointList = []
                for i in range(3, len(TopolDict[key])):
                    # If it's a list, then add that 3d point
                    pointList.append([TopolDict[key][i][0], TopolDict[key][i][1], TopolDict[key][i][2], TopolDict[key][i][3]])
                if len(pointList) < 2: numsegs = 1
                else: numsegs = 3
                ID_compartment = np.array([float(i)/numsegs for i in range(numsegs) ])*len(pointList) 
                for idc in ID_compartment:
                    xi, yi, zi, di = pointList[int(idc)]
                    h.pt3dadd(xi, yi, zi, di)
                h.pop_section()
                key.nseg = numsegs

    def set_intrinsic_properties(self, IntrinsicFactors):
        """ 
        Set intrinsic properties: set axial resistance, fill the neuron's surface with ion channels, and set their dynamical parameters

        Final factors are defined as the multiplication of the *individual* factor by the *experiment* factor, so that

            .. code-block:: python

                FinalFactors = FactorsIndividual * FactorsExperiment
                fNa, fA, fAHPs, fC, fKDR, fM, fCaL, fCaT, fHCN, fL, Ra = FinalFactors

        It is crucial to set the *individual* and *experiment* factors in the required order

        Parameters
        ----------
        IntrinsicFactors: List
            :class:`neuron_class`'s :data:`LCNhm_main.IntrinsicFactors` input: List of all intrinsic properties

            * :data:`LCNhm_configurationfile.CELLPROP_INTRINSIC_IONCHS`: Ion channels to include in the cell

            * :data:`LCNhm_configurationfile.CELLPROP_INTRINSIC_EXPERIMENT`: Additional factor multiplying the following maximum conductances and axial resistance

            * :data:`LCNhm_configurationfile.CELLPROP_INDIVIDUAL`: Number of the intrinsic genetic-algorithm `individual`

        """

        # CELLPROP_INTRINSIC_IONCHS as IonChannels
        # CELLPROP_INTRINSIC_EXPERIMENT as FactorsExperiment
        # CELLPROP_INDIVIDUAL as FactorsNum
        IonChannels, FactorsExperiment, FactorsNum = IntrinsicFactors

        # Get 'individual' factors
        with open(DirLoc+'/LCNhm-neurondata/intrinsic_individuals.txt') as file:
            FactorsIndividual = map(float, [line.split() for k, line in enumerate(file) if k==FactorsNum ][0])

        # Multiplication of the 'individual' factor by the 'experiment' factor
        FinalFactors = [ FactorsIndividual[ii] * FactorsExperiment[ii] for ii in range(len(FactorsExperiment))]

        # Separate all different factors
        fNa, fA, fAHPs, fC, fKDR, fM, fCaL, fCaT, fHCN, fL, Ra = FinalFactors

        # Somatic z-position (along radial axis)
        zSoma = h.z3d(0,sec=self.SomaList[0])
        # Intrinsic Passive Properties
        for sect in self.SectionsList:
            sect.insert('pas')
            sect.insert('extracellular')
            sect.insert('iL')
            for seg in range(sect.nseg):
                loc = float(seg+1)/(sect.nseg+1) # Location in segment
                sect.Ra = Ra # Internal resistivity (Ra), in Ohm*cm
                sect.e_pas = -50 # Membrane rest potential
                sect(loc).gl_iL = 0.0025*fL # Leak conductance
                SS = gradient_spine_scale(zDistance=h.z3d(seg,sec=sect)-zSoma, Diameter=h.diam3d(seg,sec=sect)) # Spine scale
                sect.cm = 5 if sect in self.SomaList else 1.8*5*SS  # Membrane capacitance (microFarad/cm2)
                # Specific membrane resistivity (Rm), or passive conductivity (gpas)
                if sect in self.SomaList:
                    sect(loc).g_pas = 1./gradient_membrane_resistance(0.0) # S/cm2
                elif sect in self.DendList+self.AxonList:
                    sect(loc).g_pas = SS/gradient_membrane_resistance(0.0) # S/cm2
                else:
                    sect.g_pas = SS/gradient_membrane_resistance(h.z3d(0,sec=sect)-zSoma) # S/cm2

        # HCN channels
        if 'hcn' in IonChannels:
            for sect in self.ApicList:
                sect.insert('hcn')
                sect.ehcn_hcn = -30.
                sect.Ft_hcn = 1.
                zDistance = h.z3d(0,sec=sect)-zSoma
                sect.ghcn_hcn = gradient_gHCN(zDistance)*fHCN # S/cm2
                sect.V12_hcn = gradient_V12(zDistance) # mV

        # A-Type K+ channels. Potassium conductance in S/cm2
        if 'iA' in IonChannels:
            for sect in self.SomaList+self.DendList+self.ApicList:  sect.insert('iA')
            for sect in self.SomaList: sect.gkbar_iA = 0.0025*fA
            for sect in self.DendList: sect.gkbar_iA = 0.0600*fA
            for sect in self.ApicList: sect.gkbar_iA = 0.0600*fA

        # CA2+ dependent slow AHP K+ conductance. Potassium conductance in S/cm2
        if 'iAHPs' in IonChannels:
            for sect in self.SomaList+self.DendList+self.ApicList:  sect.insert('iAHPs')
            for sect in self.SomaList: sect.gkbar_iAHPs = 0.0005*5*fAHPs
            for sect in self.DendList: sect.gkbar_iAHPs = 0.0005*fAHPs
            for sect in self.ApicList: sect.gkbar_iAHPs = 0.0005*fAHPs

        # short-duration [Ca]- and voltage-dependent current. Potassium conductance in S/cm2
        if 'iC' in IonChannels:
            for sect in self.SomaList+self.DendList+self.ApicList:  sect.insert('iC')
            for sect in self.SomaList: sect.gkbar_iC = 0.09075*fC
            for sect in self.DendList: sect.gkbar_iC = 0.03300*fC
            for sect in self.ApicList: 
                if h.z3d(0,sec=sect)-zSoma <= 200.: sect.gkbar_iC = 0.03300*fC
                elif h.z3d(0,sec=sect)-zSoma <= 350.: sect.gkbar_iC = 0.00410*fC

        # Short-duration [Ca]- and voltage-dependent current. Calcium conductance in S/cm2
        if 'iCaL' in IonChannels:
            for sect in self.SomaList:  sect.insert('iCaLs')
            for sect in self.DendList+self.ApicList:  sect.insert('iCaLd')
            for sect in self.SectionsList:  sect.insert('cad')
            for sect in self.SomaList: sect.gcalbar_iCaLs = 0.000700000*fCaL; sect.tau_cad = 50
            for sect in self.DendList: sect.gcalbar_iCaLd = 0.000031635*fCaL; sect.tau_cad = 20
            for sect in self.ApicList: sect.gcalbar_iCaLd = 0.000031635*fCaL; sect.tau_cad = 20
            for sect in self.AxonList: sect.tau_cad = 20

        # Short-duration [Ca]- and voltage-dependent current. Calcium conductance in S/cm2
        if 'iCaT' in IonChannels:
            for sect in self.SomaList+self.DendList+self.ApicList:  sect.insert('iCaT')
            for sect in self.SectionsList:  sect.insert('cad')
            for sect in self.SomaList: sect.gcatbar_iCaT = 0.00005*fCaT; sect.tau_cad = 50
            for sect in self.DendList: sect.gcatbar_iCaT = 0.00001*fCaT; sect.tau_cad = 20
            for sect in self.ApicList: sect.gcatbar_iCaT = 0.00001*fCaT; sect.tau_cad = 20
            for sect in self.AxonList: sect.tau_cad = 20

        # Delay rectifier current. Potassium conductance in S/cm2
        if 'iKDR' in IonChannels:
            for sect in self.SectionsList:  sect.insert('iKDR')
            for sect in self.SomaList: sect.gkbar_iKDR = 0.001400*fKDR
            for sect in self.AxonList: sect.gkbar_iKDR = 0.020000*fKDR
            for sect in self.DendList: sect.gkbar_iKDR = 0.000868*fKDR
            for sect in self.ApicList: sect.gkbar_iKDR = 0.000868*fKDR

        # Slowly activating voltage-dependent potassium current. Potassium conductance in S/cm2
        if 'iM' in IonChannels:
            for sect in self.SectionsList:  sect.insert('iM')
            for sect in self.SomaList: sect.gkbar_iM = 0.06*fM
            for sect in self.AxonList: sect.gkbar_iM = 0.03*fM
            for sect in self.DendList: sect.gkbar_iM = 0.06*fM
            for sect in self.ApicList: sect.gkbar_iM = 0.06*fM

        # Sodium conductance. Sodium conductance in S/cm2
        if 'iNas' in IonChannels:
            for sect in self.SectionsList:  sect.insert('iNas')
            for sect in self.SomaList: sect.gnabar_iNas = 0.007*10*fNa; sect.Frec_iNas = 1.0
            for sect in self.AxonList: sect.gnabar_iNas = 0.007*10*fNa; sect.Frec_iNas = 1.0
            for sect in self.DendList: sect.gnabar_iNas = 0.007*fNa; sect.Frec_iNas = 1.0
            for sect in self.ApicList: sect.gnabar_iNas = 0.007*fNa; sect.Frec_iNas = 1.0

    def set_synaptic_properties(self, SynapticFactors):
        """
        Set synaptic properties written in 
        *../LCHhm-neurondata/synaptic_properties.txt* (Bezaire2016)
        
        Each line of the file has a list of all the necessary
        properties to configure synapses. Let's take a line as
        an example:

            Eg.: (1) **CA3** (2) **6209** (3) **4** (4) **50** (5) **300** (6) **1.5** (7) **0** (8) **0.5** (9) **3** (10) **0.0002** (11) **0.5** (12) **276** (13) **5** (14) **3**

        These parameters are, in order:

        1. Name of input.
            Eg.: CA3 

        2. Number of NumBoutons.
            Eg.: 6209

        3. Type of ``Sections`` where to place NumBoutons.
            Eg.: 4, that is apical dendrites

        4. Minimum distance to soma (positive to apical, negative to Basal).
            Eg.: 50 um from soma, that would be proximal SR

        5. Maximum distance to soma (positive to apical, negative to Basal).
            Eg.: 300 um from soma, that would be distal SR

        6. Firing frequency (Hz).
            Eg.: 1.5 Hz

        7. Reversal potential (mV). Near 0mV would correspond to a glutamatergic input, while a -70mV to a GABAergic one.
            Eg.: 0 mV, that is excitatory

        8. Raising time constant (ms).
            Eg.: 0.5 ms

        9. Decaying time constant (ms).
            Eg.: 3.0 ms

        10. Maximum conductance (microsiemens, uS). This is later multiplied by the ``FinalFactors`` defined from the ``FactorsIndividual`` and the ``Experimental``.
             Eg.: 0.0002 uS

        11. Basal firing of the synaptic time probability distribution.
             Eg.: 0.5

        12. Theta Phase of maximum spiking probability.
             Eg.: 276 deg

        13. Left-shifting factor of the synaptic time probability distribution (A in the picture).
             Eg.: 5

        14. Right-shifting factor of the synaptic time probability distribution (B in the picture).
             Eg.: 3

        The last four parameters define the synaptic time probability distribution
        through the :func:`LCNhm_functions.synaptic_time_probability_distribution`,
        an asymmetric gaussian function. An example of the synaptic time probability
        distribution alon 10 theta cycles, and a visual definition of the
        ``Basal``, ``A`` and ``B`` parameters are shown in the image:

        .. image:: ../sphinx-docs/_images/SynapseDistribution.png

        In the case that :func:`LCNhm_configurationfile.SIMPROP_THETA_MODE` is ``False``,
        the distribution will be homogeneous along time, being (6), the firing frequency,
        the only parameter that will be taken into account.

        So for each excitatory/inhibitory input:

        1. An amount of (2) boutons are placed randomly along the neuron surface from (5) to (6)

        2. Synapse internal dynamics are set given (7), (8), (9), (10) and (11)
        
        3. Synapses are activated randomly with the probability distribution given by :func:`LCNhm_functions.synaptic_time_probability_distribution`


        Final factors are defined as the multiplication of the *individual* factor by the *experiment* factor, so that

            .. code-block:: python

                FinalFactors = FactorsIndividual * FactorsExperiment
                fNa, fA, fAHPs, fC, fKDR, fM, fCaL, fCaT, fHCN, fL, Ra = FinalFactors

        It is crucial to set the *individual* and *experiment* factors in the required order

        Parameters
        ----------
        SynapticFactors: List
            * :const:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_INPUTS` : Synaptic inputs to include in the cell

            * :const:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_EXPERIMENT` : Additional factor multiplying the following maximum conductances

            * :const:`LCNhm_configurationfile.SIMPROP_THETA_MODE` : Set theta (``True``) or not (``False``)

            * :const:`LCNhm_configurationfile.SIMPROP_THETA_PERIOD` : Theta period in milliseconds

            * :const:`LCNhm_configurationfile.SIMPROP_START_TIME` : Lapse of time before starting the simulation in milliseconds

            * :const:`LCNhm_configurationfile.SIMPROP_END_TIME` : Total duration of the simulation in milliseconds (:const:`LCNhm_configurationfile.SIMPROP_START_TIME` + :const:`LCNhm_configurationfile.SIMPROP_SIM_TIME`)

        Returns
        -------
        SynDict: Dictionary
            Dictionary with information of all synapses. For each input in :const:`LCNhm_configurationfile.CELLPROP_SYNAPTIC_INPUTS`,
            
            * ``SynSection``: NEURON ``Section`` and ``Location`` for each bouton

            * ``SynPlaces``: Three spatial coordinates for each bouton

            * ``SynTimes``: Time releases for each bouton

            is stored

        SynObjects: List
            NEURON synaptic ``Exp2Syn`` objects (go to `Exp2Syn <https://www.neuron.yale.edu/neuron/static/py_doc/modelspec/programmatic/mechanisms/mech.html?highlight=exp2syn#Exp2Syn>`_ for more information)

        NetConObjects: List
            NEURON synaptic ``NetCon`` objects (go to `NetCon <https://www.neuron.yale.edu/neuron/static/py_doc/modelspec/programmatic/network/netcon.html?highlight=netcon>`_ for more information)

        """
        # CELLPROP_SYNAPTIC_INPUTS as SynInputs
        # CELLPROP_SYNAPTIC_EXPERIMENT as FactorsExperiment
        # SIMPROP_THETA_MODE as isTheta
        # SIMPROP_THETA_PERIOD as ThetaPeriod
        # SIMPROP_START_TIME as tini
        # SIMPROP_END_TIME as tend

        # Unpack list of synaptic properties
        SynInputs, FactorsExperiment, isTheta, ThetaPeriod, tini, tend, dt = SynapticFactors
        if (tend-tini)<1.5*ThetaPeriod: warnings.warn('Warning: simulation lasts less than 1.5 times SIMPROP_THETA_PERIOD.')

        # Get 'individual' factors
        with open(DirLoc+'/LCNhm-neurondata/synaptic_individuals.txt') as file:
            FactorsIndividual = map(float, [line.split()[1:] for line in file if line.split()[0]==self.MorphoName ][0])

        # FinalFactors: Multiplication of the 'individual' factor by the 'experiment' factor
        FinalFactors = { SynInputs[ii] : (FactorsIndividual[ii]*FactorsExperiment[ii]) for ii in range(len(SynInputs)) }

        # Posible Sections into which set the synapses
        OtherApical = [ii for ii in self.ApicList if ii not in self.SomApicList]
        PosibleSections = [self.SomaList, self.AxonList, self.DendList, self.ApicList, self.SomApicList, OtherApical]

        # Position of soma
        xyzSoma = np.array([h.x3d(0,sec=self.SomaList[0]), h.y3d(0,sec=self.SomaList[0]), h.z3d(0,sec=self.SomaList[0])])
        zSoma = xyzSoma[2]

        # Initialize Returns
        SynDict = {}
        SynObjects = []
        NetConObjects = []

        if sum(FactorsIndividual)>0.:
            # Reading properties of synapses from file
            with open(DirLoc+'/LCNhm-neurondata/synaptic_properties.txt') as file:
                for line in file:

                    # Get and set parameters
                    # ----------------------
                    # Get parameters of synapses
                    SynInput = line.split()[0]
                    if SynInput in SynInputs:
                        NumBoutons, SynSection, Zmin, Zmax, FiringFreq, Erev, Tau1, Tau2, Gmax, Basal, Phase, BetaA, BetaB  = map( float, line.split()[1:] )
                        # Make new key
                        if SynInput not in SynDict.keys():
                            SynDict[SynInput] = {}
                            SynDict[SynInput]['SynSection'] = []
                            SynDict[SynInput]['SynPlaces'] = []
                            SynDict[SynInput]['SynTimes'] = []
                        # Make Int values
                        SynSection, NumBoutons = int(SynSection), int(NumBoutons)
                        # Number of synapse relases 
                        NumSynapses = FiringFreq*((tend-tini)/1000.)
                        # Multiply maximal conductance by the FinalFactor
                        Gmax *= FinalFactors[SynInput]
             
                        # Selecting places
                        # ----------------
                        # PossiblePoints: list of all [Section, Location, 3D-coords, z-coord] of the neuron between Zmin and Zmax
                        PossiblePoints = []
                        for sec in PosibleSections[SynSection-1]:
                            for seg in range(sec.nseg):
                                loc = float(seg+1)/(sec.nseg+1)
                                if (h.z3d(seg,sec=sec)-zSoma)>=Zmin and (h.z3d(seg,sec=sec)-zSoma)<=Zmax:
                                    xyz3d = np.array([h.x3d(seg,sec=sec), h.y3d(seg,sec=sec), h.z3d(seg,sec=sec)])
                                    PossiblePoints.append([sec, loc, xyz3d, xyz3d[2]])
                        PossiblePoints = np.array(PossiblePoints)
                        Zs = PossiblePoints[:,3].astype(float)
                        # Gaussian distribution around the mean ((Zmin+Zmax)/3.), and normalized to (Zmin+Zmax)/6.
                        PlaceDistribution = np.exp( -(Zs-(Zmin+Zmax)/3.)**2 / np.max([10, (Zmin+Zmax)/6.])**2)
                        PlaceDistribution /= sum(PlaceDistribution)
                        # Pick random NumBoutons PossiblePoints, according to PlaceDistribution
                        #SynPlaces = PossiblePoints[ np.in1d( Zs, np.random.choice(Zs, p=PlaceDistribution, size=NumBoutons)) ]
                        #SynPlaces = np.array(SynPlaces)
                        # Pick random NumBoutons PossiblePoints, according to the distribution, where to set the synapse (boutons)  
                        SynPlaces = [PossiblePoints[ np.where( Zs==iz )[0][0]] for iz in np.random.choice(Zs, p=PlaceDistribution, size=NumBoutons)]
                        SynPlaces = np.array(SynPlaces)

                        # Selecting times
                        # ---------------
                        # Possible times between 'tini' and 'tend'
                        PossibleTimes = np.arange(tini, tend, dt)
                        # Synaptic time distribution
                        if isTheta:
                            # If we are in theta mode, no-homogeneous distribution (given by synaptic_time_probability_distribution)
                            TimeDistribution = synaptic_time_probability_distribution(PossibleTimes, Phase, [ThetaPeriod,BetaA,BetaB])
                        else:
                            # If we are not in theta mode, homogeneous distribution
                            TimeDistribution = np.ones(len(PossibleTimes))/len(PossibleTimes)
                        # Randomness in amount of number of synapses per bouton 
                        NumSynRand = 2.*NumSynapses*np.random.random(len(SynPlaces))

                        # Setting synapses
                        # ----------------
                        for iBoutons in range(len(SynPlaces)):                        
                            # Synapse place
                            sec = SynPlaces[iBoutons][0] # Section
                            loc = SynPlaces[iBoutons][1] # Location
                            xyzSec = SynPlaces[iBoutons][2]
                            # Making GLUTAMATE/GABA A synapse at that place
                            SynObjects.append(h.Exp2Syn(sec(loc)))
                            SynObjects[-1].e = Erev # (mV)
                            SynObjects[-1].tau1 = Tau1 # (ms)
                            SynObjects[-1].tau2 = Tau2 # (ms)
                            # Taking spiking times from the above distribution
                            SynTimes = np.sort(np.random.choice(PossibleTimes, size=int(NumSynRand[iBoutons]), p=TimeDistribution))
                            # Making the synapse object
                            for tsyn in SynTimes:
                                # Connecting the synapse Object with the 'stimulator' NetCon
                                NetConObjects.append(h.NetCon(self.SomaList[0](0)._ref_v, SynObjects[-1], sec=self.SomaList[0]))
                                NetConObjects[-1].delay = tsyn # ms
                                NetConObjects[-1].threshold = -1000 # mV
                                NetConObjects[-1].weight[0] = Gmax*(1 + (np.linalg.norm(xyzSec - xyzSoma)>240)) # uS (double on distal locations)
                            # Saving for future writing
                            SynDict[SynInput]['SynSection'].append([sec,loc])
                            SynDict[SynInput]['SynPlaces'].append(xyzSec)
                            SynDict[SynInput]['SynTimes'].append(SynTimes)
                 
        return SynDict, SynObjects, NetConObjects

    def set_current_clamp(self, CurrentFactors):
        """ 
        Set current clamp

        Parameters
        ----------
        CurrentFactors: List
            :class:`neuron_class`'s :data:`LCNhm_main.IntrinsicFactors` input: List of properties for all desired current pulses
                * :const:`LCNhm_configurationfile.CURRENT_DURATION`

                List of durations of each current pulse in milliseconds

                * :const:`LCNhm_configurationfile.CURRENT_DELAY`

                List of delays of each current pulse in milliseconds

                * :const:`LCNhm_configurationfile.CURRENT_AMPLITUDES`

                List of amplitudes of each current pulse in nanoampers (nA)

                * :const:`LCNhm_configurationfile.CURRENT_SECTION`

                List of sections of each current pulse

                * :const:`LCNhm_configurationfile.CURRENT_LOCATION`

                List of location along the defined :const:`LCNhm_configurationfile.CURRENT_SECTION` of each current pulse

            Their length must be the same, the number of different current clamps

            Click any of the links for further information


        Returns
        -------
        CurrentObject: List
            :class:`neuron_class`'s :attr:`CurrentObject` attribute: List of NEURON objects containing all current clamps information
            
                [[ Clamp object 1, Duration 1, Delay 1, Amplitude 1 ],
                [  Clamp object 2, Duration 2, Delay 2, Amplitude 2 ], [...], ...... ]

            Eg.:

            .. code-block:: python

                #                   Object                 dur   del  amp 
                CurrentObject = [ [ NEURON object at soma, 100,    0,   0],
                                  [ NEURON object at soma, 200,  300,   0],
                                  ...
                                  [ NEURON object at soma, 200, 1200,   0]]
        """

        # Unpack all current parameters
        Durations, Delays, Amplitudes, Sections, Locations = CurrentFactors

        current = []
        for clamp in range(len(Durations)):
            current.append([])
            Section = Sections[clamp]
            Location = Locations[clamp]
            if   Section[:8] == 'SomaList': Section = self.SomaList[int(Section[8:])]
            elif Section[:8] == 'ApicList': Section = self.ApicList[int(Section[8:])]
            elif Section[:8] == 'DendList': Section = self.DendList[int(Section[8:])]
            elif Section[:8] == 'AxonList': Section = self.AxonList[int(Section[8:])]
            current[-1] = h.IClamp(Section(Location))
            current[-1].dur = Durations[clamp] # ms
            current[-1].delay = Delays[clamp] # ms
            current[-1].amp = Amplitudes[clamp] # nA

        return current

