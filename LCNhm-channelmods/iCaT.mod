TITLE Somatic L-type calcium channel with low threshold for activation
: used in somatic and proximal dendritic regions 
: it calculates I_Ca using channel permeability instead of conductance

COMMENT

From paper Muellner2015
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=206244&file=/CA1_multi/mechanism/cat.mod#tabs-2
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
        SUFFIX iCaT
        USEION ca READ cai,cao WRITE ica VALENCE 2
        RANGE gcatbar, minf, mtau, hinf, htau, i
}

: ----------------------------------------------------------- 

PARAMETER {                :parameters that can be entered when function is called in cell-setup 

        gcatbar = 0     (S/cm2) : initialized conductance
        ki  = 0.001     (mM)  
        cai = 5.e-5     (mM)      : initial internal Ca++ concentration
        cao = 2         (mM)      : initial external Ca++ concentration
        tfa = 1                   : time constant scaling factor
        tfi = 0.68
        tBase = 23.5    (degC)
        eca = 140       (mV)      : Ca++ reversal potential
}

: ----------------------------------------------------------- 

ASSIGNED {                       : parameters needed to solve DE

        v               (mV)
        celsius         (degC)

        gcat            (S/cm2) 
        ica             (mA/cm2)
        i               (mA/cm2)
        asdf            

        minf
        mtau            (ms)
        hinf
        htau            (ms)
        }

: ----------------------------------------------------------- 

STATE {  m h }                          : unknown parameter to be solved in the DEs 

: ----------------------------------------------------------- 

BREAKPOINT {

        SOLVE states METHOD cnexp
        gcat = gcatbar*m*m*h*h2(cai)    : maximum channel permeability
        i   = gcat*ghk(v,cai,cao)       : calcium current induced by this channel
        ica = i
}

: ----------------------------------------------------------- 

INITIAL {                               : initialize the following parameter using rates()
        rates(v)
        m = minf
        h = hinf
        gcat = gcatbar*m*m*h*h2(cai)
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

        LOCAL a

        UNITSOFF 
        a = alpm(v)
        mtau = 1.0 / (tfa*(a+betm(v)))  : estimation of activation tau
        minf = a /( a + betm(v) )       : estimation of activation steady state value

        a = alph(v)
        htau = 1.0 / (tfi*(a+beth(v)))  : estimation of inactivation tau
        hinf = a / ( a + beth( v ) )    : estimation of inactivation steady state value

        UNITSON
}

: -----------------------------------------------------------

UNITSOFF

FUNCTION h2(cai(mM)) {
        h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {

        LOCAL nu,f

        f   = KTF(celsius)/2
        nu  = v/f
        ghk = -f*(1. - (ci/co)*exp(nu))*efun(nu)
}


FUNCTION KTF(celsius (degC)) (mV) { : temperature-dependent adjustment factor
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
        if (fabs(z) < 1e-4) {
                efun = 1 - z/2
        }else{
                efun = z/(exp(z) - 1)
        }
}


FUNCTION alpm(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    alpm = 0.1967*(-1.0*v+19.88)/(exp((-1.0*v+19.88)/10.0)-1.0)
}

FUNCTION betm(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    betm = 0.046*exp(-v/22.73)
}

FUNCTION alph(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    alph = 1.6e-2*exp(-(v+57.0)/19.0)
}

FUNCTION beth(v(mV)) {
    TABLE FROM -150 TO 150 WITH 200
    beth = 1.0/(exp((-v+15)/10)+1.0)
}

UNITSON

: -----------------------------------------------------------