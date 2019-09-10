TITLE iM - slowly activating voltage-dependent potassium current

COMMENT

From M. Migliore June 2006

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

    SUFFIX iM

    USEION k READ ek WRITE ik

        : The RANGE: variables which can have a different value in each compartment
        : Variable in RANGE should also be declared in a PARAMETER or ASSIGNED block. 

    RANGE gkbar, gk, i, tau

        : GLOBAL changing its value will affect all cell

    GLOBAL minf, mtau, Q, taua, taub

        : assigned GLOBALs will be per thread

    THREADSAFE 

}

: ----------------------------------------------------------- 

: PARAMETER: Variables whose values are normally specified by the user

PARAMETER   {

    gkbar  = .06       (S/cm2)   <0,1e9>
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
    Q

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    minf
    mtau         (ms)
    taua         (ms)
    taub         (ms)
    tau          
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

        : Current iNa

    gk = gkbar*m*Q*1e-4

    i  = gk*(v - ek) 
    ik = i

    tau = mtau

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
    :if ( m < minf) { mtau = taua } else { mtau = taub }
    m' = (minf-m)/mtau

}

: ----------------------------------------------------------- 

: Functions and mathematical expressions that describe the rest of the variables

? rates
PROCEDURE rates(v(mV)) {  

    LOCAL alpha, beta
    TABLE minf, mtau DEPEND celsius, Q FROM -100 TO 100 WITH 200   


        : Computes rate and other constants at current v.
        : Call once from HOC to initialize inf at resting v.

        UNITSOFF

        :"m"

    Q = 2.3^((celsius - 23.)/10.)

    alpha = 1e-3 * vtrap(-(v+30.0), 9.0 )
    beta  = 1e-3 * vtrap( (v+30.0), 9.0 ) 

    mtau = 1.0 / ( Q * ( alpha + beta ) )
    minf = alpha / ( alpha + beta )
    
        UNITSON
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.

        if (fabs(x/y) < 1e-6) {

                vtrap = y*(1 - x/y/2)

        }else{
        
                vtrap = x/(exp(x/y) - 1)
        }
}

: ----------------------------------------------------------- 
