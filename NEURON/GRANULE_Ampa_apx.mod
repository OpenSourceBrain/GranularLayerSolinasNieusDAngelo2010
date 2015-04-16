TITLE AMPACOD  

COMMENT
	Modello AMPA dell'articolo (versione 15 settembre 2004).
ENDCOMMENT

NEURON {
	POINT_PROCESS GRANULE_Ampa_apx
	NONSPECIFIC_CURRENT i
	
	RANGE r1FIX,r2,r3,r4,r5,r1,r6,r6FIX
	RANGE g,gmax,kB,Cdur,Erev , ic
	RANGE gg1,gg2,gg3,Tdiff
	RANGE T,Tmax,Trelease 
	
	RANGE A1, A2, A3, tau_dec1, tau_dec2, tau_dec3
	RANGE tdelay,ton	 
	RANGE cvI, syntype, gmax_factor
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)	
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(um) = (micrometer)
	PI	 = (pi)	(1)
    }
    
    PARAMETER {
	syntype
	gmax_factor = 1
	: Parametri Postsinaptici
	r1FIX	= 5.4		(/ms/mM) 	:16.32	4	1				
	r2		= 0.82		(/ms)		:2.5	0.5	0.18
	r3		= 0			(/ms)		:0.1	0.18
	r4		= 0			(/ms)		:0.1	0
	r5		= 0.013		(/ms)		:0.04	0	0.063
	r6FIX	= 1.12		(/ms/mM)	:3.35	0	0
	gmax	= 400 		(nS)		: 467    1400/3   (pS)	: 3 releasing sites!
	Cdur	= 0.3		(ms)		 
	Erev	= 0			(mV)
	kB		= 0.44		(mM)
	Tmax	= 1  	(mM)	
		
	: Diffusion			
	A1 			= 0.01466	:0.0082
	A2			= 0.07126	:0.0490
	A3 			= 0.17594	:0.02558	
	tau_dec1 	= 53.6335	:80.95
	tau_dec2 	= 6.14589	:3.15
	tau_dec3 	= 1.17791	:12.56

	: Initial condition of kinetic scheme
	C_0=1
	D_0=0
	O_0=0
	cvI			= 0.1
}


ASSIGNED {
	v			(mV)		: postsynaptic voltage
	i 			(nA)		: current = g*(v - Erev)
	ic 			(nA)		: current = g*(v - Erev)
	g 			(pS)		: conductance
	r1			(/ms)
	r6			(/ms)
	T			(mM)
	Trelease	(mM)
	Tdiff		(mM)
	tdelay		(ms)
	ton			(ms)
	Tdiff_0		(mM)
	cvIsigma
	cvInoise
}

STATE {	
	C
	O
	D
	gg1
	gg2
	gg3
	sink
}	
	

INITIAL {
	C		=	C_0
	O		=	O_0
	D		=	D_0
	T		=	0 	(mM)
	Tdiff	=	0	(mM)
	Trelease=	0 	(mM)
	gg1		=	0
	gg2		=	0
	gg3		=	0   
	Tdiff_0	=	0	(mM)
	ton		=  -1   (ms)
	cvIsigma = cvI*gmax
}

PROCEDURE seed(x) { set_seed(x) }

FUNCTION SET_IC(c_0,d_0){
	C_0=c_0
	D_0=d_0
	O_0=1-C_0-D_0
	if((c_0<0)||(d_0<0)){printf("Wrong initial conditions!\n")}else{printf("Initial conditions setted to (C,D,O)=(%g,%g,%g)\n",C_0,D_0,O_0)}
}

FUNCTION SET_tdelay(R,D){ tdelay=0.25*R*R/D } :printf("Tdelay set to %g\n",tdelay)

BREAKPOINT {
	if( (t-ton)>tdelay  ){
		Tdiff=gg1+gg2+gg3
		Tdiff_0 = Tdiff
	}else{
		Tdiff=Tdiff_0+(A1+A2+A3)*(t-ton)/tdelay
	}
	Trelease	= 	T + Tdiff
	SOLVE kstates METHOD sparse
	g = (gmax+cvInoise) * O
	i = (1e-6) * g * (v - Erev) * gmax_factor
	ic = i
}


KINETIC kstates {
	: Postsynaptic scheme
	r1 = r1FIX * Trelease	:^2 / (Trelease + kB)^2
	r6 = r6FIX * Trelease	:^2 / (Trelease + kB)^2
	~ C  <-> O	(r1,r2)
	~ O  <-> D	(r3,r4)
	~ D  <-> C	(r5,r6)
	CONSERVE C+O+D = 1
	: Glutamate diffusion wave
	~ gg1 <-> sink (1/tau_dec1,0)
	~ gg2 <-> sink (1/tau_dec2,0)
	~ gg3 <-> sink (1/tau_dec3,0)
}


NET_RECEIVE(weight, on, nspike, flagtdel) {
	INITIAL {
		flagtdel=1
		nspike = 1
	}
   if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			:t0 = t
			ton=t
			on = 1		
			T=Tmax
			cvInoise=normrand(0,cvIsigma)	
		}
		net_send(Cdur, nspike)
		net_send(tdelay, flagtdel)						
    }
	if (flag == nspike) { 		
			T = 0
			on = 0
	}
	if (flag == flagtdel){
		flagtdel = flagtdel+1
		state_discontinuity(gg1,gg1+A1)	:	*x*u
		state_discontinuity(gg2,gg2+A2)	:	*x*u
		state_discontinuity(gg3,gg3+A3)	:	*x*u
	}
}	 

 
