module AvalancheModelRotatingConfigs

import QuadGK as qgk
import Roots as rts

using ..ShortDynamics


# Integral of the radiation rate over time
export int_Wncs

function int_Wncs(epsilon,weff,t,lambda=8e-7)# epsilon in a0, t in laser period, weff in omega laser
    return qgk.quadgk(time -> Wgam_t(epsilon,weff,time,lambda),0,t,rtol=1e-8)[1]
end

# Find time of emission
export solve_tem

function solve_tem(epsilon,weff,lambda=8e-7,beta=1)
    return rts.find_zero(t -> int_Wncs(epsilon,weff,t,lambda)-beta,(1e-5,10),xatol=1e-3)
end

# Parameters of the charge at emission time
export chi_em

function chi_em(epsilon,weff,lambda=8e-7,beta=1)#Quantum parameter at tem(beta)
    return chi_e(epsilon,weff,solve_tem(epsilon,weff,lambda,beta),lambda)
end


export gamma_em

function gamma_em(epsilon,weff,lambda=8e-7,beta=1)#Lorentz factor at tem(beta)
    return gamma_e(epsilon,solve_tem(epsilon,weff,lambda,beta))
end

# Estimation of the growth rate
export GR_steady_state

function GR_steady_state(Wrad,Wcr,Nu) #formula of the growth rate from quasi steady state assumption, rates in 1/tau
    delta = (Wcr-Nu)*(Wcr-Nu) - 4*Wcr*(Nu - 2*Wrad)
    return 0.5*(sqrt(delta)-Wcr -Nu)
end


export GR_CPSW

function GR_CPSW(epsilon,weff,Nu,lambda=8e-7,beta=1)#epsilon in a0, weff in omega nu in 1/tau
    tem = solve_tem(epsilon,weff,lambda,beta)
    Wrad = Wgam_t(epsilon,weff,tem,lambda)
    Wcr = Wpair_t(epsilon,weff,tem,lambda)
    return GR_steady_state_mig(Wrad,Wcr,Nu)
end


export GR_pure_rotating

function GR_pure_rotating(epsilon,lambda=8e-7,beta=1)# Pure rotating field is a particular case of CPSW with nu = 0 and weff = 0.5omega
    return GR_CPSW(epsilon,0.5,0,lambda,beta)
end


export GR_CPSW_eps2weff

function GR_CPSW_eps2weff(epsilon,eps2weff,Nu,lambda=8e-7,beta=1)#epsilon in a0, weff in omega nu in 1/tau
    weff = eps2weff/(epsilon*epsilon)
    tem = solve_tem(epsilon,weff,lambda,beta)
    Wrad = Wgam_t(epsilon,weff,tem,lambda)
    Wcr = Wpair_t(epsilon,weff,tem,lambda)
    return GR_steady_state_mig(Wrad,Wcr,Nu)
end

end