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

**Note**: the package is not yet registered in the official julia package repositories, which is why the github url has to be specified to install it.

## Documentation

The package is mainly used to be able to predict avalanches growth rates in strong electromagnetic field configurations. The functions used to do these estimate are presented in the following. After the main section, all of the intermediate functions used to achieve the growth rates computations are described.

### Growth Rates predictions

### Strong-Field QED processes rates

Here are documented the functions needed for SFQED rates computations.

#### Nonlinear Breit-Wheeler pair production

* `Wnbw(gamma,chi,lambda)` returns the nonlinear Breit-Wheeler rate value for a gamma photon with normalized energy `gamma` and quantum parameter `chi`. The result is normalized to the corresponding frequency of a wave with wavelength `lambda` : $c/\lambda$. If the wavelength `lambda` is not specified, the value $\lambda = 0.8 \,  \mu\rm{m}$ is used.

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

* `Wncs(gamma,chi,lambda)` returns the nonlinear Compton scattering rate value for a lepton  with Lorentz factor `gamma` and quantum parameter `chi`. The result is normalized to the corresponding frequency of a wave with wavelength `lambda` : $c/\lambda$. If the wavelength `lambda` is not specified, the value $\lambda = 0.8 \,  \mu\rm{m}$ is used.

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

### Invariants and effective frequency
