TITLE NMDA sinaptico

COMMENT
	Modello NMDA dell articolo (versione 15 settembre 2004).
	C'e' solo il Tdiff e non T (pulse)
ENDCOMMENT

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GRANULE_Nmda_apx
	NONSPECIFIC_CURRENT i

	RANGE Rb,Ru,Rd,Rr,Ro, Rc,rb
	RANGE g,gmax,Cdur,Erev, ic 
			
	RANGE MgBlock,v0_block,k_block
	RANGE gg1,gg2,gg3,Tdiff
	RANGE T,Tmax,Trelease 
	RANGE A1, A2, A3, tau_dec1, tau_dec2, tau_dec3
	RANGE tdelay,ton
	RANGE cvI, syntype, gmax_factor
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(uM) = (micro/liter)
	PI	= (pi)		(1)
    }
    
    PARAMETER {
	syntype
	gmax_factor = 1
	Rb			=  5		(/ms/mM)  	:15: binding  
	Ru			=  1		(/ms)		:15e-3	: unbinding
	Rd			=  12e-3  	(/ms)		: desensitization
	Rr			=  9e-4		(/ms)		: resensitization 
	Ro			=  3e-2 	(/ms)		: opening
	Rc			=  0.966	(/ms)		: closing
	Erev		= -3.7  	(mV)	
	gmax		= 5333.3	(pS)	
	v0_block 	= -20 		(mV)	
	k_block 	= 13		(mV)
	Tmax		= 0.5  		(mM)
	Cdur		= 0.3		(ms)

	: Diffusion			
	A1 			= 0.01466	:0.0082
	A2			= 0.07126	:0.0490
	A3 			= 0.17594	:0.02558	
	tau_dec1 	= 53.6335	:80.95
	tau_dec2 	= 6.14589	:3.15
	tau_dec3 	= 1.17791	:12.56
						
	: Initial condition of kinetic scheme
	C0_0=1
	C1_0=0
	C2_0=0
	D_0=0
	O_0=0
	
	cvI			= 0.1
	celsius (degC)
}


ASSIGNED {
	v			(mV)		: postsynaptic voltage
	i 			(nA)		: current = g*(v - Erev)
	ic 			(nA)		: current = g*(v - Erev)
	g 			(pS)		: actual conductance
	rb			(/ms)    : binding
	MgBlock
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
	C0		: unbound
	C1		: single bound
	C2		: double bound
	D		: desensitized
	O		: open
	gg1
	gg2
	gg3
	sink
}

INITIAL {
	rates(v)
	C0		=	C0_0
	C1		=	C1_0
	C2		=	C2_0
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

FUNCTION SET_tdelay(R,D){ tdelay=0.25*R*R/D } :printf("Tdelay set to %g\n",tdelay)

BREAKPOINT {
	rates(v)
: 	Tdiff		=	gg1+gg2+gg3	
: 	Trelease	=	T + Tdiff
	
	if( (t-ton)>tdelay) {
		Tdiff=gg1+gg2+gg3
		Tdiff_0 = Tdiff
	}else{
		Tdiff=Tdiff_0+(A1+A2+A3)*(t-ton)/tdelay
	}
	Trelease	= 	Tdiff :T + Tdiff
	SOLVE kstates METHOD sparse
	g = (gmax+cvInoise) * O 			
	i = (1e-6) * g * (v - Erev) * MgBlock * gmax_factor
	ic = i
}

KINETIC kstates {	
	rb = Rb * Trelease 
	~ C0 <-> C1	(rb,Ru) 	: (fattore*rb,Ru) qui 2* per descrizione part.identiche
	~ C1 <-> C2	(rb,Ru)		: (rb,fattore*Ru)	idem
	~ C2 <-> D	(Rd,Rr)
	~ C2 <-> O	(Ro,Rc)
	CONSERVE C0+C1+C2+D+O = 1
	: Glutamate diffusion wave
	~ gg1 <-> sink (1/tau_dec1,0)
	~ gg2 <-> sink (1/tau_dec2,0)
	~ gg3 <-> sink (1/tau_dec3,0)
}

PROCEDURE rates(v(mV)) {
	: E' necessario includere DEPEND v0_block,k_block per aggiornare le tabelle!
	TABLE MgBlock DEPEND v0_block,k_block FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + exp ( - ( v - v0_block ) / k_block ) )
}


NET_RECEIVE(weight, on, nspike, flagtdel) {
	INITIAL {
		flagtdel=1
		nspike = 1
	}
   if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			ton = t
			on = 1		
			:T=Tmax
			cvInoise=normrand(0,cvIsigma)
			:printf("%g\n",cvInoise)
		}
		net_send(Cdur, nspike)
		net_send(tdelay, flagtdel)						
    }
	if (flag == nspike) { 
			:t0 = t
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
