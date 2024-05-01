# PairAvalanchesQED.jl

PairAvalanchesQED.jl is a [Julia](https://julialang.org/) package providing numerical tools to study electron-positron pairs avalanches in the context of strong-field quantum electrodynamics and ultra intense electromagnetic fields.

The present code was created and used during the work published [here](https://arxiv.org/abs/2402.04225). It was gathered in the form of a Julia package in order to be shared easily. Please do not hesitate to share it, use it in your own work or research, or build your own new tools from it, as it is shared under a GPL-3.0 license. If you do so, please make sure to cite the publication of origin above. The following documentation will use conventions and notations from the above publication, it is highly recommended to take a look at it to better understand the documentation and how to use the package.

## How to install

To install this package, first make sure to have [Julia](https://julialang.org/) 1.9 or above on your computer.

Then in a Julia [interactive session](https://docs.julialang.org/en/v1/manual/getting-started/) (REPL) execute the two following lines

```julia
using Pkg
Pkg.add(url="https://github.com/amercurib/PairAvalanchesQED.jl")
```

After that the PairAvalanchesQED package can be called like any other package like this

```julia
using PairAvalanchesQED
```

or like this

```julia
import PairAvalanchesQED as avalanches
```

**Note**: the package is not yet registered in the official julia  repositories, which is why the github url has to be specified to install it.

## Documentation

The package is mainly used to be able to predict avalanches growth rates in strong electromagnetic field configurations. The functions used to do these estimate are presented in the following. After the main section, all of the intermediate functions used to achieve the growth rates computations are described. In the following $\omega$ is the reference (laser) frequency and $\tau$ the (laser) period such that $\omega = 2 \pi c/\tau$. 

### Table of contents
1. [Growth Rates predictions](https://github.com/amercurib/PairAvalanchesQED.jl#growth-rates-predictions)
2. [Strong-Field QED processes rates](https://github.com/amercurib/PairAvalanchesQED.jl#strong-field-qed-processes-rates)
    1. [Nonlinear Breit-Wheeler pair production](https://github.com/amercurib/PairAvalanchesQED.jl#nonlinear-breit-wheeler-pair-production)
    2. [Nonlinear inverse Compton scattering](https://github.com/amercurib/PairAvalanchesQED.jl#nonlinear-inverse-compton-scattering)
3. [Short time dynamics model](https://github.com/amercurib/PairAvalanchesQED.jl#short-time-dynamics-model)
4. [Invariants and effective frequency](https://github.com/amercurib/PairAvalanchesQED.jl#invariants-and-effective-frequency)

### Growth Rates predictions 
The main functions of this modules are the following, which can be used to compute growth rate estimates. the default values of the functions for optional parameters are set to be the same as in the publication.


* `GR_CPSW(epsilon,weff,Nu=1,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)`
  returns the growth rate prediction for a circularly polarised standing wave normalized to $1/\tau$. 
    * `epsilon` is the invariant $\epsilon$, should be given with the same normalization as the laser field strength $a_0$.
    * `weff` is the effective frequency $\omega_{\rm eff}$ in units of $\omega$.
    * `Nu` is the migration frequency in units of $1/\tau$ default value is `1`.
    * `lambda` the reference laser wavelength, by default $\lambda = 0.8 \,  \mu\rm{m}$.
    * `chi_switch`is the value of quantum parameter at emission $\chi_{\textrm{em}}$ above which the results is given from the high field model, the default value is `10`.
    * `beta` the $\beta$ constant in the equation of $t_{\textrm{em}}$, by default `1`.
    * `tem_precision` the precision on characteristic time of emission $t_{em}$, by default `1e-4`.

* `GR_pure_rotating(epsilon,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)`
  returns the growth rate prediction for a purely rotating electric field, normalized to $1/\tau$. The parameters are the same as for `GR_CPSW` except for `epsilon` and `weff` which are not necessary here.

* `GR_CPSW_esp2weff(epsilon,esp2weff,Nu=1,lambda=8e-7,chi_switch=10,beta=1,tem_precision=1e-4)`
  does the same as `GR_CPSW` but using $\epsilon^2 \omega_{\rm eff}$ instead of $\omega_{\rm eff}$.

The following are intermediate or more specific functions for growth rate estimates.

* `GR_CPSW_HighIntensity(epsilon,weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the growth rate for a circularly polarized standing wave estimate using the high intensity part of the model, the result is normalized de $1/\tau$. The parameters are the same as for `GR_CPSW` except for `chi_switch` which is not needed here.

* `GR_CPSW_MidIntensity(epsilon,weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the growth rate for a circularly polarized standing wave estimate using the lower and intermediate intensity part of the model, the result is normalized de $1/\tau$. The parameters are the same as for `GR_CPSW` except for `chi_switch` which is not needed here.

* `GR_CPSW_HighIntensity_eps2weff(epsilon,eps2weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)`
  does the same as `GR_HighIntesity` but using $\epsilon^2 \omega_{\rm eff}$ instead of $\omega_{\rm eff}$. 

* `GR_CPSW_MidIntensity_eps2weff(epsilon,eps2weff,Nu=1,lambda=8e-7,beta=1,tem_precision=1e-4)`
  does the same as `GR_midIntesity` but using $\epsilon^2 \omega_{\rm eff}$ instead of $\omega_{\rm eff}$.

* `GR_pure_rotating_HighIntensity(epsilon,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the growth rate prediction for a purely rotating electric field normalized to $1/\tau$, using the high intensity part of the model. The parameters are the same as for `GR_pure_rotating` except for `chi_switch`which is not needed here.

* `GR_pure_rotating_MidIntensity(epsilon,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the growth rate prediction for a purely rotating electric field normalized to $1/\tau$, using the lower and intermediate intensity part of the model. The parameters are the same as for `GR_pure_rotating` except for `chi_switch`which is not needed here.

* `GR_steady_state(Wrad,Wcr,Nu)`
  returns the value of the growth rate using the formula obtained from the quasi-steady state equations $\Gamma = \dfrac{W_{\mathrm{cr}} + \nu}{2}\left[\sqrt{1 + \frac{4W_{\mathrm{cr}}(2W_{\mathrm{rad}} - \nu)}{(W_{\mathrm{cr}} + \nu)^2}} -1 \right]$.

* `solve_tem(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the characteristic time of emission $t_{em}$ as given by the equation $\int_0^{t_{\textrm{em}}}  W_{\textrm{CS}}(t') \textrm{d} t' = \beta$.
    * `epsilon` is the invariant $\epsilon$, should be given with the same normalization as the laser field strength $a_0$.
    * `weff` is the effective frequency $\omega_{\rm eff}$ in units of $\omega$
    * `lambda` the reference laser wavelength, by default $\lambda = 0.8 \,  \mu\rm{m}$.
    * `beta`the $\beta$ constant in the equation of $t_{\textrm{em}}$, by default `1`.
    * `tem_precision` the precision on $t_{em}$, by default `1e-4`.

* `chi_em(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the quantum parameter for a lepton initially at rest at the characteristic time of emission $t_{em}$. The parameters are the same as for `solve_tem`.

* `gamma_em(epsilon,weff,lambda=8e-7,beta=1,tem_precision=1e-4)`
  returns the Lorentz factor for a lepton initially at rest at the characteristic time of emission $t_{em}$. The parameters are the same as for `solve_tem`.

* `int_Wncs(epsilon,weff,t,lambda=8e-7)`
  returns the integral of the nonlinear Compton scattering rate over time `t` for a lepton, using classical trajectories given by the short time dynamics model. The invariant `epsilon` should be normalized as the laser field strength $a_0 = \dfrac{e E}{mc\omega}$ and the effective frequency `weff` should be in units of $\omega$. If no wavelength `lambda` is specified the function will use a reference wavelength of $0.8 \mu \textrm{m}$.

  

### Strong-Field QED processes rates 

Here are documented the functions needed for SFQED rates computations.

#### Nonlinear Breit-Wheeler pair production

* `Wnbw(gamma,chi,lambda)` returns the nonlinear Breit-Wheeler rate value for a gamma photon with normalized energy `gamma` and quantum parameter `chi`. The result is normalized to the frequency $1/\tau$ of a wave with wavelength `lambda` : $c/\lambda$. If the wavelength `lambda` is not specified, the value $\lambda = 0.8 \,  \mu\rm{m}$ is used.

* `Wnbw_SI(gamma,chi)`
returns the nonlinear Breit-Wheeler rate value in SI units ($\textrm{s}^{-1}$) for a gamma photon with normalized energy `gamma` and quantum parameter `chi`.

* `b0(chi)`
returns the value of the function $b_0(\chi_{\gamma}) = \gamma_{\gamma} W_{\rm BW}(\chi_{\gamma},\gamma_{\gamma})/W_0$ needed for the computation of the $W_{\rm BW}$ rate. The function has been tabulated and if the quantum parameter `chi` is between `0.01`and `10000` the function will return the value interpolated from the table. Outside of this range, it will return the values from the asymptotic formulas given by the functions `b0_high_chi(chi)` and `b0_low_chi(chi)`.

* `b0_low_chi(chi)`
returns the value of the $\chi_{\gamma} \ll 1$ asymptote given by $b_0(\chi_{\gamma}) \simeq 0.344 \times \chi_{\gamma} \exp(-\frac{8}{3\chi_{\gamma}})$.

* `b0_high_chi(chi)`
returns the value of the $\chi_{\gamma} \gg 1$ asymptote given by $b_0(\chi_{\gamma}) \simeq 0.569 \times \chi_{\gamma}^{2/3}$.

* `b0_from_integral(chi)`
returns the value of the function $b_0(\chi_{\gamma})$ without using the tables and asymptotes like `b0(chi)`, but by computing the integral. This function is used to generates the table `b0(chi)`, and is much slower than the version using table.

* `compute_b0_table(filename)`
    generate a table for `b0`as a `.jld` file at the directory specified by `filename`. Note that the values of quantum parameter in the table are logarithmically spaced.

    By default the parameters of this table are as used by the package. However more arguments can be passed to tune the output: `compute_b0_table(filename,Npoints_table,min_chi,max_chi,period_show_progress)`
  * `Npoints_table`: number of points of the table, by default `30000`.
  * `min_chi`: minimum quantum parameter of the table, by default `0.01`.
  * `max_chi`: maximum quantum parameter of the table, by default `10000`.
  * `period_show_progress`: number of table points after which to print the computation status, by default `1000`.

#### Nonlinear inverse Compton scattering

* `Wncs(gamma,chi,lambda)` returns the nonlinear Compton scattering rate value for a lepton  with Lorentz factor `gamma` and quantum parameter `chi`. The result is normalized to the frequency $1/\tau$ of a wave with wavelength `lambda` : $c/\lambda$. If the wavelength `lambda` is not specified, the value $\lambda = 0.8 \,  \mu\rm{m}$ is used.

* `Wncs_SI(gamma,chi)`
returns the nonlinear Compton scattering rate value in SI units ($\textrm{s}^{-1}$) for a lepton with Lorentz factor `gamma` and quantum parameter `chi`.

* `c0(chi)`
returns the value of the function $c_0(\chi_{e}) = \gamma_{e} W_{\rm CS}(\chi_{e},\gamma_{e})/W_0$ needed for the computation of the $W_{\rm CS}$ rate. The function has been tabulated and if the quantum parameter `chi` is between `0.001`and `10000` the function will return the value interpolated from the table. Outside of this range, it will return the values from the asymptotic formulas given by the functions `c0_high_chi(chi)` and `c0_low_chi(chi)`.

* `c0_low_chi(chi)`
returns the value of the $\chi_{e} \ll 1$ asymptote given by $c_0(\chi_{e}) \simeq 2.16 \times \chi_{e} $.

* `c0_high_chi(chi)`
returns the value of the $\chi_{e} \gg 1$ asymptote given by $c_0(\chi_{e}) \simeq 2.19 \times \chi_{e}^{2/3}$.

* `c0_from_integral(chi)`
returns the value of the function $c_0(\chi_{e})$ without using the tables and asymptotes like `c0(chi)`, but by computing the integral. This function is used to generates the table `c0(chi)`, and is much slower than the version using table.

* `compute_c0_table(filename)`
    generate a table for `c0`as a `.jld` file at the directory specified by `filename`. Note that the values of quantum parameter in the table are logarithmically spaced.

    By default the parameters of this table are as used by the package. However more arguments can be passed to tune the output: `compute_c0_table(filename,Npoints_table,min_chi,max_chi,period_show_progress)`
  * `Npoints_table`: number of points of the table, by default `35000`.
  * `min_chi`: minimum quantum parameter of the table, by default `0.001`.
  * `max_chi`: maximum quantum parameter of the table, by default `10000`.
  * `period_show_progress`: number of table points between each print of the computation status, by default `1000`.

### Short time dynamics model

* `chi_e(epsilon, weff, t)`
  returns the value of the quantum parameter of an electron or positron at time `t`, following the approximation at short times $ \chi_e(t) \simeq \dfrac{\epsilon^2 \omega_{\rm eff}}{E_S^2\tau_C} t^2$. The particle is considered initially at rest. The parameter `epsilon` should be given with the same normalization as the laser field strength $a_0$, `weff` should be given in units of $\omega$ and `t`in laser periods $\tau$.

* `gamma_e(epsilon, t)`
  returns the value of the Lorentz factor of an electron or positron at time `t`, following the approximation at short times $ \gamma_{e}(t)  \simeq \dfrac{e\epsilon  t}{mc}$. The particle is considered initially at rest. The parameter `epsilon` should be given with the same normalization as the laser field strength $a_0$ and `t`in laser periods $\tau$.

* `chi_e_SI(epsilon, weff, t)`
  does the same as `chi_e ` but parameters should be given in SI units.
  
* `gamma_e_SI(epsilon, t)`
  does the same as `\gamma_e` but parameters should be given in SI units.

* `Wgam_t(epsilon,weff,t)`
  returns the rate of the gamma photons emission by nonlinear Compton scattering by a lepton with the Lorentz factor and the quantum parameter given by the short time dynamics model at time `t`. The rate is given in units of $1/\tau$ and `epsilon`and `weff`should have the same units as in `chi_e`.

* `Wpair_t(epsilon,weff,t)`
  returns the rate of the pair creation by nonlinear Breit-Wheeler for a photon, with the normalized energy and the quantum parameter being the same as the one of a lepton in the short time dynamics model at time `t`. The rate is given in units of $1/\tau$ and `epsilon` and `weff` should have the same units as in `chi_e`.

* `Wgam_t_SI(epsilon,weff,t)`
  does the same as `Wgam_t` but the parameters and the results are in SI units.

* `Wpair_t_SI(epsilon,weff,t)`
  does the same as `Wpair_t` but the parameters and the results are in SI units.



### Invariants and effective frequency

This part of the package contains the functions to compute the $\omega_{\textrm{eff}}$ and $\epsilon$ for a given field configuration, in order to use them in the short time dynamics model and then in the growth rate estimates. The main function of use in this model is the following

* `weff_epsilon_from_field(EB_T,dx,dy,dz)`
  returns $\omega_{\textrm{eff}}$, $\epsilon$ and $\epsilon^{2}\omega_{\textrm{eff}}$ values over space. The parameter `EB_T` should be a 4 dimensional array, the first three dimensions being the coordinates in space, and the last one the components of the electric and magnetic field (in this order) at the corresponding space coordinate. The parameters `dx`, `dy` and `dz` are the resolutions over the `EB_T` array along the three directions.

In the following will be explained the intermediate functions used to build `weff_epsilon_from_field`.

* `computing_invariants(EB_T)`
  returns the Lorentz invariants $\mathcal{F}$ and $\epsilon$ over space for a given field configuration. `EB_T`is the same as for `weff_epsilon_from_field`.

* `make_TF(EB_T)`
  returns the components electromagnetic field strength tensor $F^{\mu}_{\phantom{1} \nu}$ corresponding to the field configuration `EB_T`. The results is therefore a 5 dimensional array, where the first three indices are the space coordinates and the last two are the component of the field strength tensor `EB_T` is the same as for `weff_epsilon_from_field`.

* `make_square_field_tensor(TF)`
  returns the square matrix of the components of the electromagnetic tensor, given the tensor components over space `TF`.

* `computing_derivatives_EB(EB_T,dx,dy,dz)`
  returns the spatial derivatives and time derivatives (using Maxwell's equations) of the components of the field over space contained in `EB_T`. The results has therefore on more dimension than `EB_T`where the indices corresponds to the time and space derivatives (in this order). The parameters are the same as for `weff_epsilon_from_field`.

* `make_derivative_tensor(DEB_T)`
  given the derivatives of the fields over space `DEB_T` from the `computing_derivatives_EB` function, build the derivatives of the electromagnetic tensor components $F^{\mu}_{\phantom{1} \nu , \sigma}$.

* `make_J_matrix(Eps_inv,TF2)`
  * computes the $J$ matrix over space, and returns its values in a 5 dimensional array. `Eps_inv` is a 4 dimensional array of $\epsilon$ over space and `TF2` is the square matrix of the electromagnetic field components as given by the function `make_square_field_tensor`.

* `invert_J_matrix(J,F_inv)`
  invert the $J$ matrix for all points of space where `F_inv`is strictly positive.

* `computing_eigenvector_f1(Eps_inv,EB_T)`
  compute the eigen 4-vector components corresponding to the eigen value $\epsilon$ over space. `Eps_inv` is given by `computing_invariants` and `EB_T`is the same as for `weff_epsilon_from_field`.

* `lower_f1_index(f1)`
  lower the index of the `f1` four vector over space.


* `compute_weff(F_inv,DTF,Jm1,f1_low,f1)`
  returns the value of $\omega_{\textrm{eff}}$ over space given the intermediate quantities provided by the previously described functions.
    * `F_inv` is the $\mathcal{F} invariant over space,
    * `DTF` is the array of the derivative of the electromagnetic tensor components over space,
    * `Jm1` is the inverse of the $J$ matrix over space,
    * `f1_low` is the array of the lowered components of the `f1` eigen 4-vector over space,
    * `f1` is the array of the components of the `f1` eigen 4-vector over space.


