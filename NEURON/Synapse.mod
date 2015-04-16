TITLE SynCO

COMMENT
Solinas et al. 2007, Frontiers in Cellular Neuroscience.
Synaptic model C=O gating scheme.
tau_1, tau_2 are the binding and unbinding time constants.
ENDCOMMENT

NEURON {
	POINT_PROCESS Synapse
	NONSPECIFIC_CURRENT i	
	RANGE g,gmax,kB,Cdur,Erev 
	RANGE tau_1,tau_2
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)	
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(um) = (micrometer)
	PI   = (pi)(1)
}

PARAMETER {
	: Parametri Postsinaptici
	tau1		= .4		(/ms/mM) 	 
        tau2		= 3		(/ms)
        gmax		= 1150 		(pS)
	Cdur		= 0.3		(ms)		 
	Erev		= 0		(mV)
	kB		= 0.44		(mM)		

         : Parametri Presinaptici
	tau_1 		= 3 (ms) 	< 1e-9, 1e9 >
	tau_rec 	= 35.1 (ms) 	< 1e-9, 1e9 > 	
	tau_facil 	= 10.8 (ms) 	< 0, 1e9 > 	

	U 		= 0.1 (1) 	< 0, 1 >
	u0 		= 0 (1) 	< 0, 1 >	: se u0=0 al primo colpo y=U
	Tmax		= 1  (mM)
}


ASSIGNED {
	v		(mV)		: postsynaptic voltage
	i 		(nA)		: current = g*(v - Erev)
	g 		(pS)		: conductance
	r1		(/ms)
	r6		(/ms)
	T		(mM)
        Trelease        (mM)
	ton		(ms)	
	
	x
	PRE
}

STATE {	
	C
	O
}	
	

INITIAL {
	C		=	1
	O		=	0
	T		=	0 	(mM)
	Trelease	=	0 	(mM)
	ton		=  	-1   (ms)
	PRE		=	0
}

BREAKPOINT {
	Trelease=T
	SOLVE kstates METHOD sparse
	g =gmax * O
	i = (1e-6) * g * (v - Erev) : 1e-6 fA=>nA
}


KINETIC kstates {
	: Postsynaptic scheme
	r1 = 1/tau1 * T
	~ C  <-> O	(r1,1/tau2)
	CONSERVE C+O = 1
}


NET_RECEIVE(weight, on, nspike, t0 (ms),y, z, u, tsyn (ms)) {
	INITIAL {
		nspike = 1
		y=0
		z=0
		u=u0
		tsyn=t
	}
   	if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			ton=t
			t0=t
			on=1				
			z=z*exp(-(t-tsyn)/tau_rec)
			z=z+(y*(exp(-(t - tsyn)/tau_1)-exp(-(t-tsyn)/tau_rec))/(tau_1/tau_rec-1))
			y=y*exp(-(t-tsyn)/tau_1)			
			x=1-y-z
			if(tau_facil>0){ 
				u=u*exp(-(t-tsyn)/tau_facil)
				u=u+U*(1-u)							
			}else{u=U}
			y=y+x*u
			PRE=y
			   T=Tmax :*y
			tsyn=t						
		}
		net_send(Cdur,nspike)
    	}
	if(flag==nspike){ 		
		T = 0
		on = 0
	}
}	 

 
