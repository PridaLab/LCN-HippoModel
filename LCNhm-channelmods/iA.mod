TITLE iA - A-type K+ clannels
 
COMMENT

From paper Cutsuridis2015

ENDCOMMENT

: ----------------------------------------------------------- 
 
UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
}

: ----------------------------------------------------------- 

: NEURON: to insert different values in each compartment.

NEURON  {
    
        : Name of the distributed mechanism 

    SUFFIX iA

    USEION k READ ek WRITE ik

        : The RANGE: variables which can have a different value in each compartment
        : Variable in RANGE should also be declared in a PARAMETER or ASSIGNED block. 

    RANGE gkbar, gk, i

        : GLOBAL changing its value will affect all cell

    GLOBAL minf, mtau, hinf

        : assigned GLOBALs will be per thread

    THREADSAFE 

}

: ----------------------------------------------------------- 

: PARAMETER: Variables whose values are normally specified by the user

PARAMETER   {

    gkbar  = .0025      (S/cm2)   <0,1e9>
    ek     = -80        (mV)
    mtau    = 0.1       (ms)
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
    hinf
    htau         (ms)
}

: ----------------------------------------------------------- 

: STATE: dependent variables, or unknowns in differential equations, families of algebraic equations, or kinetic reaction schemes
: These variables do not need to be declared in the NEURON block.   

STATE { 
    
    m h

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

    gk = gkbar*m*h
    i  = gk*(v - ek) 
    ik = i

}

: ----------------------------------------------------------- 

: INITIAL: initialization of all STATEs

INITIAL     {

    rates(v)
    m = minf
    h = hinf

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
    h' = (hinf-h)/htau

}

: ----------------------------------------------------------- 

: Functions and mathematical expressions that describe the rest of the variables

? rates
PROCEDURE rates( v(mV) ) {  

    LOCAL alpha, beta
    TABLE minf, hinf, htau FROM -100 TO 100 WITH 200   

        : Computes rate and other constants at current v.
        : Call once from HOC to initialize inf at resting v.

        UNITSOFF

        :"m"
    
    alpha = 0.01 * vtrap( -(v+21.3) , 35.0 )
    beta  = 0.01 * vtrap(  (v+21.3) , 35.0 )
    minf  = alpha/(alpha+beta)

        :"h"
    
    alpha = -0.01 * vtrap( (v+58.0) , 8.2 )
    beta  = -0.01 * vtrap(-(v+58.0) , 8.2 )
    hinf  = alpha/(alpha+beta)

    htau = htauv(v)

        UNITSON

}

: ----------------------------------------------------------- 

FUNCTION htauv( v (mV) ) {  

        UNITSOFF

        : Traps for 0 in denominator of rate eqns.
    
    if ( v > - 20.0 ) {

            htauv = 5. + 2.6*( v+20. )/10.

    }else{

            htauv = 5.
    
    }

        UNITSON
}

FUNCTION vtrap(x,y) {  :Traps for 0 in denominator of rate eqns.

        if (fabs(x/y) < 1e-6) {

                vtrap = y*(1.0 - x/y/2.0)

        }else{
        
                vtrap = x/(exp(x/y) - 1.0)
        }
}