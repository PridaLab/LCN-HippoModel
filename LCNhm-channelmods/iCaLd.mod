TITLE Dendritic L-type calcium channel with low threshold for activation
: used in somatic and proximal dendritic regions 
: it calculates I_Ca using channel permeability instead of conductance

COMMENT

From paper Muellner2015
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=206244&file=/CA1_multi/mechanism/calH.mod#tabs-2

ENDCOMMENT

UNITS {
        (mA)    = (milliamp)
        (mV)    = (millivolt)
        (molar) = (1/liter)
        (mM)    = (millimolar)
        (S)     = (siemens)
        FARADAY = 96520 (coul)
        R       = (k-mole) (joule/degC)     : (k-mole) is stored as the product of Boltzmanns constant and Avogadros number. 
        KTOMV   = .0853 (mV/degC)
}

: ----------------------------------------------------------- 

NEURON {
        SUFFIX iCaLd
        USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf, mtau, hinf, htau, i
}

: ----------------------------------------------------------- 

PARAMETER {                :parameters that can be entered when function is called in cell-setup 

        gcalbar = 0     (S/cm2)   : initialized conductance
        ki  = 0.001     (mM)  
        cai = 5.e-5     (mM)      : initial internal Ca++ concentration
        cao = 2         (mM)      : initial external Ca++ concentration
        tfa = 5                   : time constant scaling factor
        eca = 140       (mV)      : Ca++ reversal potential
}

: ----------------------------------------------------------- 

ASSIGNED {                       : parameters needed to solve DE

        v               (mV)
        celsius         (degC)

        gcal            (S/cm2) 
        ica             (mA/cm2)
        i               (mA/cm2)

        minf
        mtau            (ms)
        hinf
        htau            (ms)
        }

: ----------------------------------------------------------- 

STATE { 
        m
        h 
}                      : unknown parameter to be solved in the DEs 

: ----------------------------------------------------------- 

BREAKPOINT {

        SOLVE states METHOD cnexp
        gcal = gcalbar*m*m*m*h
        i    = gcal*(v - eca)
        ica  = i
}

: ----------------------------------------------------------- 

INITIAL {                        : initialize the following parameter using rates()

        rates(v)
        m   = 0
        h   = 1
}

: ----------------------------------------------------------- 

? states
DERIVATIVE states {
    
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau

}

: ----------------------------------------------------------- 

PROCEDURE rates(v (mV)) { :callable from hoc

        UNITSOFF 

        minf = 1.0 / (1.0 + exp((v+37.0)/(-1.0)))         : Ca activation 
        mtau = 3.6

        hinf = 1.0 / (1.0 + exp((v+41.0)/(0.5)))          : Ca inactivation
        htau = 29.0

        UNITSON
}

: -----------------------------------------------------------
