
module Rates

### Physical constants for computations (SI units)
e = 1.602176620e-19 #elementary charge
hbar = 1.054571628e-34 #reduced Planck constant
c = 299792458 #speed of light

alpha = 7.2973525664e-3 #fine structure constant
me = 9.109382616e-31 #electron mass
E_mass_e = me*c*c #electron rest energy
W0 = 2*alpha*me*c*c/(3*hbar) #constant rate multiplicative factor appearing in SF processes rates
TW0 = 1/W0 #corresponding time of W0 rate
re = 2.8179e-15 #classical radius of the electron


import SpecialFunctions as spe
import QuadGK as qgk
import JLD as jld

### Nonlinear Breit-Wheeler process

b0_high_factor= 45*3^(2/3)*spe.gamma(2/3)^4/(56*pi*pi) #constant for high chi asymptot
b0_low_factor = 3/16*(3/2)^(1.5) #constant for low chi asymptot

#import the table to speed up the computation of b0
path_table_b0 = joinpath(pkgdir("PairAvalanchesQED"),"tables","b0_table.jld")
Table_b0 = jld.load(path_table_b0,"b0")
chi_table_b0 = jld.load(path_table_b0,"chi")
size_table_b0 = length(Table_b0)
chi_min_b0 = chi_table_b0[1]
chi_max_b0 = last(chi_table_b0)
log_chi_min_b0 = log10(chi_min_b0)
log_chi_max_b0 = log10(chi_max_b0)


export b0_low_chi

function b0_low_chi(chi)
    return b0_low_factor*chi*exp(-8/(3*chi))
end


export b0_high_chi

function b0_high_chi(chi)
    return b0_high_factor*chi^(2/3)
end


export b0

function b0(chi)
    if chi >= chi_max_b0
        return b0_high_chi(chi)
    elseif chi <= chi_min_b0 
        return b0_low_chi(chi)
    else
        i = (size_table_b0-1)*(log10(chi) - log_chi_min_b0)/(log_chi_max_b0 -log_chi_min_b0) +1 #+1 because arrays start at index 1 in Julia
        idown = floor(Int64,i)
        iup = idown+1
        return (Table_b0[iup]-Table_b0[idown])/(chi_table_b0[iup]-chi_table_b0[idown])*(chi-chi_table_b0[idown]) + Table_b0[idown] 
    end
end


export Wnbw_SI

function Wnbw_SI(gamma,chi)
    return W0*b0(chi)/gamma
end


export Wnbw

function Wnbw(gamma,chi,lambda=8e-7) # gives the result in terms of laser frequency 1/tau, lambda is wavelength
    return W0*b0(chi)/gamma*lambda/c
end



#To compute tables 



function f_pair(xi,chi) # This integrated in xi will give b0
    mu = 2/(3*chi*xi*(1-xi))
    return sqrt(3)/(2*pi*xi*(1-xi))*( 2*(1-2*xi)*spe.besselk(5/3,mu)/(3*chi*(1-xi)) + spe.besselk(2/3,mu)  )
end


export b0_from_integral

function b0_from_integral(chi) #integrated f_rad for xi from 0 to 1
    return qgk.quadgk(xi -> f_pair(xi,chi), 0, 1, rtol=1e-8)[1]
end


####!!!!! Put function to make the table here !!!!!####




### Nonlinear Compton Scattering

c0_low_factor = sqrt(3)/(4*pi)*3*spe.gamma(1/6)*spe.gamma(11/6) #constant for low chi asymptot
c0_high_factor = sqrt(3)/(4*pi)*3^(2/3)*spe.gamma(2/3)*(2*spe.gamma(5/3)*spe.gamma(1/3) + spe.gamma(2/3)*spe.gamma(7/3)/2) #constant for high chi asymptot

#import the table to speed up the computation of c0
path_table_c0 = joinpath(pkgdir("PairAvalanchesQED"),"tables","c0_table.jld")
Table_c0 = jld.load(path_table_c0,"c0")
chi_table_c0 = jld.load(path_table_c0,"chi")
size_table_c0 = length(Table_c0)
chi_min_c0 = chi_table_c0[1]
chi_max_c0 = last(chi_table_c0)
log_chi_min_c0 = log10(chi_min_c0)
log_chi_max_c0 = log10(chi_max_c0)

export c0_low_chi

function c0_low_chi(chi)
    return c0_low_factor*chi
end


export c0_high_chi

function c0_high_chi(chi)
    return c0_high_factor*chi^(2/3)
end


export c0

function c0(chi)
    if chi >= chi_max_c0
        return c0_high_chi(chi)
    elseif chi <= chi_min_c0 
        return c0_low_chi(chi)
    else
        i = (size_table_c0-1)*(log10(chi) - log_chi_min_c0)/(log_chi_max_c0 -log_chi_min_c0) +1 #+1 because arrays start at index 1 in Julia
        idown = floor(Int64,i)
        iup = idown+1
        return (Table_c0[iup]-Table_c0[idown])/(chi_table_c0[iup]-chi_table_c0[idown])*(chi-chi_table_c0[idown]) + Table_c0[idown] 
    end
end


export Wncs_SI

function Wncs_SI(gamma,chi) # returns nonlinear Compton emission rate in SI units
    return W0*c0(chi)/gamma
end


export Wncs

function Wncs(gamma,chi,lambda=8e-7) # gives the result in terms of laser frequency 1/tau, lambda is wavelength
    return W0*c0(chi)/gamma*lambda/c
end

#To compute tables

function f_rad(xi,chi) # This integrated in xi will give c0
    mu = 2*xi/(3*chi*(1-xi))
    return sqrt(3)/(2*pi)*(mu*spe.besselk(5/3,mu)/(1-xi) + xi*xi*spe.besselk(2/3,mu)/(1-xi) )
end


export c0_from_integral

function c0_from_integral(chi) #integrated f_rad for xi from 0 to 1
    return qgk.quadgk(xi -> f_rad(xi,chi), 0, 1, rtol=1e-8)[1]
end

####!!!!! Put function to make the table here !!!!!####

end