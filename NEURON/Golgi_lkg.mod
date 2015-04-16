TITLE Cerebellum Golgi Cell Model

COMMENT

Author: S. Solinas
Data from: Solinas et al. Frontiers in Cellular Neuroscience 2007
Last revised: April 2006

ENDCOMMENT

NEURON {
	SUFFIX Golgi_lkg
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,gbar_Q10
	RANGE el, gbar, g	, ic
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v (mV)
	Q10_diff	= 1.5
	gbar = 21e-6 (mho/cm2)
	celsius  (degC)
	el = -55 (mV)
}

ASSIGNED { 
	i (mA/cm2) 
	gbar_Q10 (mho/cm2)
	ic
	g
}

BREAKPOINT { 
	gbar_Q10 = gbar*(Q10_diff^((celsius-23)/10))
	i = gbar_Q10 * (v - el ) 
	ic = i
	g = gbar_Q10
}
