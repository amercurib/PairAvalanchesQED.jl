module ShortDynamics

using ..Rates # import rates functions

### Physical constants for computations (SI units)
e = 1.602176620e-19 #elementary charge
hbar = 1.054571628e-34 #reduced Planck constant
c = 299792458 #speed of light
me = 9.109382616e-31 #electron mass





#gamma and chi estimates from short time dynamics model


export chi_e_SI

function chi_e_SI(epsilon,weff,t)
    return e*e*hbar/(me^3*c^4)*epsilon*epsilon*weff*t*t
end


export gamma_e_SI

function gamma_e_SI(epsilon,t)
    return e*epsilon*t/(me*c)
end


export chi_e

function chi_e(epsilon,weff,t,lambda=8e-7)#time should be given in laser period, epsilon in a0, omega_eff in omega
    omega = 2*pi*c/lambda
    return 4*pi^2*hbar*epsilon*epsilon*weff*omega*t*t/(me*c^2)
end


export gamma_e

function gamma_e(epsilon,t) # time should be given in laser period  , epsilon in a0
    return 2*pi*epsilon*t
end


#Rates in function of time computed with gamma and chi estimates


export Wgam_t_SI

function Wgam_t_SI(epsilon,weff,t)
    return Wncs(gamma_e_SI(epsilon,t),chi_e_SI(epsilon,weff,t))
end


export Wpair_t_SI

function Wpair_t_SI(epsilon,weff,t)
    return Wnbw(gamma_e_SI(epsilon,t),chi_e_SI(epsilon,weff,t))
end


export Wgam_t

function Wgam_t(epsilon,weff,t,lambda=8e-7) #time should be given in laser period, epsilon in a0, omega_eff in omega
    return Wncs(gamma_e(epsilon,t),chi_e(epsilon,weff,t,lambda),lambda) 
end


export Wpair_t

function Wpair_t(epsilon,weff,t,lambda=8e-7) #time should be given in laser period, epsilon in a0, omega_eff in omega
    return Wnbw(gamma_e(epsilon,t),chi_e(epsilon,weff,t,lambda),lambda) 
end

end