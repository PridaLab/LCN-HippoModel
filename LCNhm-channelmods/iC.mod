TITLE iC - short-duration [Ca]- and voltage-dependent current
 
COMMENT

From paper Cutsuridis2015

ENDCOMMENT

: ----------------------------------------------------------- 
 
UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
    (molar) = (1/liter)
    (mM)    = (millimolar)
    R       = (k-mole) (joule/degC)     : (k-mole) is stored as the product of Boltzmanns constant and Avogadros number. 
    FARADAY = (faraday) (coulombs)      : (faraday) is stored in Coul/mole
}

: ----------------------------------------------------------- 

: NEURON: to insert different values in each compartment.

NEURON  {
    
        : Name of the distributed mechanism 

    SUFFIX iC

    USEION k READ ek WRITE ik
    USEION ca READ cai

        : The RANGE: variables which can have a different value in each compartment
        : Variable in RANGE should also be declared in a PARAMETER or ASSIGNED block. 

    RANGE gkbar, gk, i

        : GLOBAL changing its value will affect all cell

    GLOBAL minf, mtau

        : assigned GLOBALs will be per thread

    THREADSAFE 

}

: ----------------------------------------------------------- 

: PARAMETER: Variables whose values are normally specified by the user

PARAMETER   {

    gkbar = .09075   (S/cm2)   <0,1e9>
    ek     = -80        (mV)
}

: ----------------------------------------------------------- 

: ASSIGNED: declare two kinds of variables: 
:   - those that are given values outside the mod file
:   - those that appear on the left hand side within the mod file
: They will not be visible at the hoc level unless it is declared in RANGE or GLOBAL

ASSIGNED {
    
    v            (mV)
    celsius      (degC)

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    cai          (mM)

    minf
    mtau         (ms)
}

: ----------------------------------------------------------- 

: STATE: dependent variables, or unknowns in differential equations, families of algebraic equations, or kinetic reaction schemes
: These variables do not need to be declared in the NEURON block.   

STATE { 
    
    m 

}

: ----------------------------------------------------------- 

: BREAKPOINT: Main computation block of the mechanism. 
: Remark: by the end of the BREAKPOINT block, all variables are consistent with the new time. 
: If a mechanism has STATEs, this block must contain one SOLVE statement that tells how the values of the STATEs will be computed over each time step. 

? currents
BREAKPOINT  {

        : Calculate variables of the STATE block 

    SOLVE states METHOD cnexp

        : Current ik

    gk = gkbar*m

    i  = gk*(v - ek)
    ik = i 

}

: ----------------------------------------------------------- 

: INITIAL: initialization of all STATEs

INITIAL     {

    rates(v,cai)
    m = minf

}

: ----------------------------------------------------------- 

: Compute derivatives of the STATEs that are described by differential equations (y' = expression)
: Equations are integrated using the numerical method specified by the SOLVE statement in the BREAKPOINT block. 
:   - cnexp: is appropriate when expression=f(y,x) is linear in y and involves no other states 
:   - rate(v): assigns values to the voltage sensitive parameters of this equation.

? states
DERIVATIVE states   {
    
    rates(v,cai)
    m' = (minf-m)/mtau

}

: ----------------------------------------------------------- 

: Functions and mathematical expressions that describe the rest of the variables

? rates
PROCEDURE rates( v(mV) , cai (mM) ) {  

    LOCAL alpha, beta, Q

        UNITSOFF

        : Computes rate and other constants at current v.
        : Call once from HOC to initialize inf at resting v.

    Q = FARADAY / ( R * ( 273.16 + celsius ) )

    alpha = ( 0.48 (/ms) ) / ( 1 + ( 0.18 (mM) / cai )*exp( -1.68*v*Q ) )
    beta  = ( 0.28 (/ms) ) / ( 1 + cai/( 0.011 (mM) * exp( -2*v*Q ) ) )

        :"m"
    
    mtau = 1/(alpha+beta)
    minf = alpha/mtau

        UNITSON

}

: ----------------------------------------------------------- 
