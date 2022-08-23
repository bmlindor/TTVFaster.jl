"""
    TTVFaster

Computes first order eccentricity transit timing variations (TTVs) with respect to the following initial properties: the planet-star mass ratio [μ], the initial transit time (of the averaged orbit) [t0], the mean orbital period [Per], the eccentricity [e], and the longitude of periastron [barred ω].
"""

module TTVFaster

include("ttv_wrapper.jl")
export Planet_plane_hk, ttv_wrapper, compute_ttv!
# export Planet_plane

end
