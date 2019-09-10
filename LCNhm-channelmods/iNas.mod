TITLE iNas.mod  - sodium channels
 
COMMENT
 This is the original Hodgkin-Huxley treatment for the set of sodium, 
  potassium, and leakage channels found in the squid giant axon membrane.
  ("A quantitative description of membrane current and its application 
  conduction and excitation in nerve" J.Physiol. (Lond.) 117:500-544 (1952).)
 Membrane voltage is in absolute mV and has been reversed in polarity
  from the original HH convention and shifted to reflect a resting potential
  of -65 mV.
 Remember to set celsius=6.3 (or whatever) in your HOC file.
 See squid.hoc for an example of a simulation using this model.
 SW Jaslove  6 March, 1992
ENDCOMMENT
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
	(S) = (siemens)
    R       = (k-mole) (joule/degC)     : (k-mole) is stored as the product of Boltzmanns constant and Avogadros number. 
    FARADAY = (faraday) (coulombs)      : (faraday) is stored in Coul/mole
}
 
? interface
NEURON {
    SUFFIX iNas
    USEION na READ ena WRITE ina
    RANGE gnabar, gna, Frec, i
    GLOBAL minf, hinf, sinf, mtau, htau, stau
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
    gnabar = .12 (S/cm2)	<0,1e9>
    Frec = 0.
}
 
STATE {
    m h s
}
 
ASSIGNED {
    v       (mV)
    celsius (degC)
    ena     (mV)

	gna     (S/cm2)
    ina     (mA/cm2)
    i       (mA/cm2)
    minf 
    hinf 
    sinf
	mtau    (ms) 
    htau    (ms) 
    stau    (ms)
}
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar*m*m*m*h*s
	i   = gna*(v - ena)
    ina = i
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
    s = sinf
}

? states
DERIVATIVE states {  
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    s' = (sinf-s)/stau
}
 

? rates
PROCEDURE rates(v(mV)) {  :Computes rate and other constants at current v.
                  :Call once from HOC to initialize inf at resting v.

    LOCAL  alpha, beta, Q
    TABLE minf, mtau, hinf, htau, sinf, stau DEPEND celsius FROM -100 TO 100 WITH 200

        UNITSOFF

    : from Magee1995, Hoffman1997, modified

            :"m" sodium activation system
    alpha = 0.182*vtrap( -(v+45), 4.5 )
    beta = 0.124*vtrap( (v+45), 4.5 )	
    mtau = 0.8/(alpha+beta)
    minf = alpha/(alpha+beta)

            :"h" sodium inactivation system
    alpha = 0.08*vtrap(-(v+40),3.0)
    beta  = 0.0005*vtrap( (v+10.0), 5. )
	htau = 1/(alpha+beta)
    hinf = 1.0 / ( 1 + exp((v+58)/5.) )

            :"s" sodium slow variable
    Q = FARADAY / ( R * ( 273.16 + celsius ) )
    sinf = ( 1. + Frec*exp( (v+58.)/2. ) )/( 1. + exp( (v+58.)/2. ) )
    stau = 0.00333 (ms) * exp( 0.0024 (/mV) * (v+60.0) * Q )/( 1.0 + exp( 0.0012 (/mV) * (v+60.0) * Q )  )
}
 
FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.
    if (fabs(x/y) < 1e-6) {
            vtrap = y*(1 - x/y/2)
    }else{
            vtrap = x/(exp(x/y) - 1)
    }
}
 
        UNITSON
