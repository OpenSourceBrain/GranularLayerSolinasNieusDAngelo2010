TITLE AMPACOD  

COMMENT
	Modello AMPA dell articolo (versione 15 settembre 2004).
	Adattato per approx il deterministico
ENDCOMMENT

NEURON {
	POINT_PROCESS AMPA_PF_GO_nodiff
	NONSPECIFIC_CURRENT i	
	RANGE r1FIX,r2,r3,r4,r5,r1,r6,r6FIX
	RANGE g,gmax,kB,Cdur,Erev 
	RANGE T,Tmax,Trelease 	
	RANGE tau_1,tau_rec,tau_facil,U,u0	
	RANGE tdelay,ton,syntype	 
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
    syntype = 0
    : Parametri Postsinaptici
	r1FIX		= 5.4		(/ms/mM) 	 
        r2		= 0.8		(/ms)		: 0.82
	r3		= 0		(/ms)		 
	r4		= 0		(/ms)		 
        r5		= 0.014		(/ms)		: 0.013
       r6FIX		= 1.2		(/ms/mM)	: 1.12
        gmax		= 2000 		(nS)		: 1200
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
         Trelease          (mM)
	ton		(ms)	
	
	x
	PRE
}

STATE {	
	C
	O
	D
}	
	

INITIAL {
	C		=	1
	O		=	0
	D		=	0
	T		=	0 	(mM)
	Trelease	=	0 	(mM)
	ton		=  	-1   (ms)
	PRE		=	0
}

BREAKPOINT {
	Trelease=Tmax*U : T
	SOLVE kstates METHOD sparse
	g =gmax * O
	i = (1e-6) * g * (v - Erev) 
}


KINETIC kstates {
	: Postsynaptic scheme
	r1 = r1FIX * Trelease^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease^2 / (Trelease + kB)^2
	~ C  <-> O	(r1,r2)
	~ O  <-> D	(r3,r4)
	~ D  <-> C	(r5,r6)
	CONSERVE C+O+D = 1
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
			   T=Tmax :*y this removes the presyn dynamics
			tsyn=t						
		}
		net_send(Cdur,nspike)
    	}
	if(flag==nspike){ 		
		T = 0
		on = 0
	}
}	 

 
