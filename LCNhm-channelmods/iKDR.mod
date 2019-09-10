TITLE iKDR - delay rectifier current
 
COMMENT

From paper Cutsuridis2015

ENDCOMMENT

: ----------------------------------------------------------- 
 
UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
    R       = (k-mole) (joule/degC)     : (k-mole) is stored as the product of Boltzmanns constant and Avogadros number. 
    FARADAY = (faraday) (coulombs)      : (faraday) is stored in Coul/mole
}

: ----------------------------------------------------------- 

: NEURON: to insert different values in each compartment.

NEURON  {
    
        : Name of the distributed mechanism 

    SUFFIX iKDR

    USEION k READ ek WRITE ik

        : The RANGE: variables which can have a different value in each compartment
        : Variable in RANGE should also be declared in a PARAMETER or ASSIGNED block. 

    RANGE gkbar, gk, i

        : GLOBAL changing its value will affect all cell

    GLOBAL minf

        : assigned GLOBALs will be per thread

    THREADSAFE 

}

: ----------------------------------------------------------- 

: PARAMETER: Variables whose values are normally specified by the user

PARAMETER   {

    gkbar  = .0014      (S/cm2)   <0,1e9>
    ek     = -80        (mV)
    mtau   = 3.50      (ms)
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

    minf
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

    gk = gkbar*m*m
    i  = gk*(v - ek) 
    ik = i

}

: ----------------------------------------------------------- 

: INITIAL: initialization of all STATEs

INITIAL     {

    rates(v)
    m = minf

}

: ----------------------------------------------------------- 

: Compute derivatives of the STATEs that are described by differential equations (y' = expression)
: Equations are integrated using the numerical method specified by the SOLVE statement in the BREAKPOINT block. 
:   - cnexp: is appropriate when expression=f(y,x) is linear in y and involves no other states 
:   - rate(v): assigns values to the voltage sensitive parameters of this equation.

? states
DERIVATIVE states   {
    
    rates(v)
    m' = (minf-m)/mtau

}

: ----------------------------------------------------------- 

: Functions and mathematical expressions that describe the rest of the variables

? rates
PROCEDURE rates(v(mV)) {  

    TABLE minf FROM -100 TO 100 WITH 200   

        : Computes rate and other constants at current v.
        : Call once from HOC to initialize inf at resting v.

        UNITSOFF

        :"m"
    
    minf = 1.0/( 1.0 + exp( -(v+42.0)/3.0 ) )

        UNITSON

}

: ----------------------------------------------------------- 
