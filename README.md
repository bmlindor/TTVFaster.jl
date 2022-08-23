# TTVFaster
First order eccentricity transit timing variations (TTVs) computed in Agol &amp; Deck (2015)
[![ADS](https://img.shields.io/badge/ADS-2016APJ...818...177A-blue)](https://ui.adsabs.harvard.edu/abs/2016ApJ...818..177A/abstract) [![arXiv](https://img.shields.io/badge/arXiv-1509.01623-brightgreen)](http://arxiv.org/abs/1509.01623)

This implements equation (33) from that paper by computing the Laplace
coefficients using a series solution due to Jack Wisdom, computing
the f_{1,j}^{(+-k)} coefficients given in equation (34) using the functions u and
v_+- with coefficients given in Table 1.

## Installation
You can install the registered TTVFaster repo as a Julia package with the `Pkg` manager.
- the repo from the package registry has been tested on Julia v1.3.0  
```julia
julia> using Pkg Pkg.add("TTVFaster.jl")
```
In its current state, the package computes the TTVs of a multi-transiting planetary system where _two_ planets are observed to be transiting. 
If you intend to modify the source code for more than 2 transiting planets, please create a GitHub fork to develop your own version. 
- make sure to replace `your-GitHub-username` with your actual GitHub username in the code below
```julia
julia> Pkg.develop(PackageSpec(url="git@github.com:your-GitHub-username/TTVFaster.jl.git"))
```
## Usage
TTVFaster computes TTVs with respect to 5 properties for each planet: \mu,t0,Period, e cos(omega), e sin(omega);
 where \mu is the mass ratio of the planet to the star, t0 is the initial transit time (of the averaged orbit), 
Period is the mean orbital period, e is the eccentricity, and omega is the longitude of periastron. 

### Example
The file kepler62ef_planets.txt in the examples/ directory contains
a comma-separated set of 10 parameters that describe a system with two planets similar to Kepler-62e/f. 

``` julia
julia> using TTVFaster,DelimitedFiles

julia> data=readdlm("kepler62ef_planets.txt",',',Float64)  
1x10 Array{Float64,2}:
 3.02306e-5  122.386  -16.5926  -0.00127324  0.0026446  1.67874e-5  267.307  155.466  -0.0025544  0.00117917

julia> include("test_ttv.jl")  
test_ttv (generic function with 4 methods)

julia> @time ttv1,ttv2=test_ttv(5,40,20,data);  # inputs are jmax,ntrans1,ntrans2,data
  0.982326 seconds (2.04 M allocations: 98.466 MiB, 12.08% gc time)
julia> @time ttv1,ttv2=test_ttv(5,40,20,data);  
  0.001171 seconds (331 allocations: 21.922 KiB)
```

This computes the TTVs and writes them to the files inner_ttv.txt and outer_ttv.txt. 
Note that the TTVs are stored in the variables ttv1 and ttv2, as well. 
The test_ttv.jl routine accepts jmax (the maximum j to sum to, in this example 5),
ntrans1 (number of transits of the inner planet), ntrans2 (the number of transits of the outer planet), 
and data which contains the parameters of both planets.
