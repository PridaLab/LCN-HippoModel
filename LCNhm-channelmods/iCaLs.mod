TITLE Somatic L-type calcium channel with low threshold for activation
: used in somatic and proximal dendritic regions 
: it calculates I_Ca using channel permeability instead of conductance

COMMENT

From paper Muellner2015
https://senselab.med.yale.edu/modeldb/ShowModel.cshtml?model=206244&file=/CA1_multi/mechanism/cal.mod#tabs-2

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
        SUFFIX iCaLs
        USEION ca READ cai,cao WRITE ica
        RANGE gcalbar, minf, mtau, i
}

: ----------------------------------------------------------- 

PARAMETER {                :parameters that can be entered when function is called in cell-setup 

        gcalbar = 0     (S/cm2) : initialized conductance
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
        }

: ----------------------------------------------------------- 

STATE {  m  }                      : unknown parameter to be solved in the DEs 

: ----------------------------------------------------------- 

BREAKPOINT {

        SOLVE states METHOD cnexp
        gcal = gcalbar*m*h2(cai) : maximum channel permeability
        i   = gcal*ghk(v,cai,cao): calcium current induced by this channel
        ica = i
}

: ----------------------------------------------------------- 

INITIAL {                        : initialize the following parameter using rates()
        rates(v)
        m = minf
        gcal = gcalbar*m*h2(cai)
}

: ----------------------------------------------------------- 

? states
DERIVATIVE states {
    
    rates(v)
    m' = (minf-m)/mtau

}

: ----------------------------------------------------------- 

PROCEDURE rates(v (mV)) { :callable from hoc

        LOCAL a

        UNITSOFF 
        a = alpm(v)
        mtau = 1.0 / (tfa*(a+betm(v)))  : estimation of activation tau
        minf = a / (a+betm(v))          : estimation of activation steady state value

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
        alpm = - 0.055 * (-27.01+v )/(exp( -(-27.01+v)/3.8) - 1)
}


FUNCTION betm(v(mV)) {
        TABLE FROM -150 TO 150 WITH 200
        betm =0.94*exp(-(63.01+v)/17)
}

UNITSON

: -----------------------------------------------------------