module AvalancheModelRotatingConfigs

import QuadGK as qgk
import Roots as rts

using ..ShortDynamics
using ..Rates


### Physical constants for computations (SI units)
e = 1.602176620e-19 #elementary charge
hbar = 1.054571628e-34 #reduced Planck constant
c = 299792458 #speed of light

me = 9.109382616e-31 #electron mass
ES = me^2. * c^3. / (e*hbar) #Schwinger field


# Integral of the radiation rate over time
export int_Wncs

function int_Wncs(epsilon,weff,t,lambda=8e-7)# epsilon in a0, t in laser period, weff in omega laser
    return qgk.quadgk(time -> Wgam_t(epsilon,weff,time,lambda),0,t,rtol=1e-8)[1]
end

# Find time of emission
export solve_tem

function solve_tem(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)
    return rts.find_zero(t -> int_Wncs(epsilon,weff,t,lambda)-beta,(1e-5,10),xatol=tem_precision)
end

# Parameters of the charge at emission time
export chi_em

function chi_em(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)#Quantum parameter at tem(beta)
    return chi_e(epsilon,weff,solve_tem(epsilon,weff,lambda,beta,tem_precision),lambda)
end


export gamma_em

function gamma_em(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)#Lorentz factor at tem(beta)
    return gamma_e(epsilon,solve_tem(epsilon,weff,lambda,beta,tem_precision))
end

# Estimation of the growth rate
export GR_steady_state

function GR_steady_state(Wrad,Wcr,Nu) #formula of the growth rate from quasi steady state assumption, rates in 1/tau
    delta = (Wcr-Nu)*(Wcr-Nu) - 4*Wcr*(Nu - 2*Wrad)
    return 0.5*(sqrt(delta)-Wcr -Nu)
end


# Estimate at mid and low intensities

export GR_CPSW_MidIntensity

function GR_CPSW_MidIntensity(epsilon,weff,Nu,lambda=8e-7,beta=1,tem_precision=1e-4)#epsilon in a0, weff in omega nu in 1/tau
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    Wrad = Wgam_t(epsilon,weff,tem,lambda)
    Wcr = Wpair_t(epsilon,weff,tem,lambda)
    return GR_steady_state_mig(Wrad,Wcr,Nu)
end


export GR_pure_rotating_MidIntensity

function GR_pure_rotating_MidIntensity(epsilon,lambda=8e-7,beta=1,tem_precision=1e-4)# Pure rotating field is a particular case of CPSW with nu = 0 and weff = 0.5omega
    return GR_CPSW(epsilon,0.5,0,lambda,beta,tem_precision)
end


export GR_CPSW_eps2weff_MidIntensity

function GR_CPSW_eps2weff(epsilon,eps2weff,Nu,lambda=8e-7,beta=1,tem_precision=1e-4)#epsilon in a0, weff in omega nu in 1/tau
    weff = eps2weff/(epsilon*epsilon)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    Wrad = Wgam_t(epsilon,weff,tem,lambda)
    Wcr = Wpair_t(epsilon,weff,tem,lambda)
    return GR_steady_state_mig(Wrad,Wcr,Nu)
end


#Estimate a high intensity


export GR_CPSW_HighIntensity

function GR_CPSW_HighIntensity(epsilon,weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)
    omega = 2*pi*c/lambda
    a0S = e *ES/(me*c*omega)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    Wgam = Wgam_t_normT(epsilon,weff,tem,lambda)
    Wpair = epsilon/a0S * Wnbw_normT(1,1,lambda)
    return GR_steady_state_mig(Wgam,Wpair,Nu)
end


export GR_pure_rotating_HighIntensity

function GR_pure_rotating_HighIntensity(epsilon,lambda=8e-7,beta=1,tem_precision=1e-4)
    GR_CPSW_HighIntensity(epsilon,0.5,0,lambda,beta,tem_precision)
end


export GR_pure_rotating_HighIntensity_eps2weff

function GR_CPSW_HighIntensity_eps2weff(epsilon,eps2weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)
    weff = eps2weff/(epsilon*epsilon)
    omega = 2*pi*c/lambda
    a0S = e *ES/(me*c*omega)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    Wgam = Wgam_t_normT(epsilon,weff,tem,lambda)
    Wpair = epsilon/a0S * Wnbw_normT(1,1,lambda)
    return GR_steady_state_mig(Wgam,Wpair,Nu)
end


# Overall model combining all intensities


export GR_CPSW

function GR_CPSW(epsilon,weff,Nu=1,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    chi_em = chi_e(epsilon,weff,tem,lambda)
    if chi_em < chi_switch
        return GR_CPSW_MidIntensity(epsilon,weff,Nu,lambda,beta,tem_precision)
    else 
        return GR_CPSW_HighIntensity(epsilon,weff,Nu,lambda,beta,tem_precision)
    end
end


export GR_pure_rotating

function GR_pure_rotating(epsilon,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    chi_em = chi_e(epsilon,weff,tem,lambda)
    if chi_em < chi_switch
        return GR_pure_rotating_MidIntensity(epsilon,lambda,beta,tem_precision)
    else 
        return GR_pure_rotating_HighIntensity(epsilon,lambda,beta,tem_precision)
    end
end


export GR_CPSW_eps2weff

function GR_CPSW_esp2weff(epsilon,esp2weff,Nu=1,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)
    tem = solve_tem(epsilon,weff,lambda,beta,tem_precision)
    chi_em = chi_e(epsilon,weff,tem,lambda)
    if chi_em < chi_switch
        return GR_CPSW_MidIntensity_esp2weff(epsilon,esp2weff,Nu,lambda,beta,tem_precision)
    else 
        return GR_CPSW_HighIntensity_esp2weff(epsilon,esp2weff,Nu,lambda,beta,tem_precision)
    end
end



end