NEURON {
	SUFFIX BK
	USEION k READ ek WRITE ik VALENCE 1
	USEION cal2 READ cal2i VALENCE 2
	RANGE gbar, minf, mtau
        GLOBAL shift
        GLOBAL tau_factor

        RANGE i_output, output
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(uM)	= (micromolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
        gbar = 1e-6 (S/cm2)
        shift = 0 (mV)
        tau_factor = 1 (ms)
        q10    = 2.5
        slope = 11.1 (/mV)
        mtau_min = 0.01
        
}

ASSIGNED {
        v       (mV)

        ek (mV)
        ik (mA/cm2)

        i_output
        output


        
		minf
		mtau 	(ms)
		cal2i	(mM)
		celsius	(degC)
		vhalf	(mV)

        q
}

STATE {
        m 
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        output    = m*gbar
        i_output  = output*(v-ek)
	ik = i_output
}
 
INITIAL {
		rates(v, cal2i)
		m = minf
                output    = m*gbar
                i_output  = output*(v-ek)
		ik = i_output
                q = q10 ^ ((celsius - 23) / 10)

}

DERIVATIVE states {  
        rates(v, cal2i)
        m' = (minf-m)/mtau
}


FUNCTION alpha(v, ca) { LOCAL alpha_max, k, vh
    alpha_max = 0 + 5.46548878 / (1 + exp(-(log10(ca) + 1.61322084) / 0.24917329)) 
    k = 8.263928565595934 + 1.2552204 / (1 + exp(-(log10(ca) + 1.61382309) / 0.03180212))  
    vh = 52.094347461863926 - 69.4428835 / (1 + exp(-(log10(ca) + 1.70836237) / 0.0236872607))
    alpha = alpha_max / (1 + exp(-(v - vh) / k))
}

FUNCTION beta(v, ca) { LOCAL vh, k
    vh = -51.449688490633456 + 0.17726174 * exp(-(log10(ca) - 4.07832681) / 1.02895756)
    k = -24.653642743189497
    
    if (fabs((v - vh)/k) < 1e-5) {
        beta = 1.0
    } else {
        beta = -1.0 / ( exp(-(v - vh) / k) - 1) * (v - vhalf) / k
    }

}

FUNCTION func_vhalf(cai) { LOCAL ca_conc_log
  ca_conc_log = log10(cai) + 3
  if (ca_conc_log <= -1) {
    func_vhalf = 152.0
  } else if (ca_conc_log >= 3) {
    func_vhalf = -47.7
  } else {
     func_vhalf = -49.925 * ca_conc_log + 102.075
  }
}


PROCEDURE rates(v(mV), cai (mM)) { LOCAL a, b, vsh
              vsh = v + shift
                                   
              a  = alpha(vsh, cai)
              b  = beta(vsh, cai)
                                   
              : definition of tau from the literature                              
              mtau = tau_factor / (a + b)

              if(mtau < mtau_min) {
                mtau = mtau_min
              }
                                   
              mtau = mtau / q
                                   
              minf = 1 / (1 + exp(-(vsh - func_vhalf(cai))/slope))
}


 
UNITSON 





