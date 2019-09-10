TITLE iL.mod - leak channels
 
UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
	(S) = (siemens)
}
 
? interface
NEURON {
    SUFFIX iL
    NONSPECIFIC_CURRENT il
    RANGE gl, el
	THREADSAFE : assigned GLOBALs will be per thread
}
 
PARAMETER {
    gl = .0025 (S/cm2)	<0,1e9>
    el = -54.3 (mV)
}
 
ASSIGNED {
    v (mV)
    celsius (degC)

    il (mA/cm2)
}
 
? currents
BREAKPOINT {

    il = gl*(v - el)
}
