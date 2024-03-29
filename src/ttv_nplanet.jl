include("compute_ttv.jl")
"""
ttv_nplanet(nplanet,jmax,ntrans,params)

Computes TTVs with TTVFaster for N planets with pairwise TTV calculation.

Arguments:
  nplanet: Number of planets
  jmax:    Maximum j over which to sum the TTV calculation
  ntrans:  Number of transits for each planet
  params:  Parameters of each planet 
Returns:
  ttvs:    Computed transit timing variations for the transiting planets.
"""

function ttv_nplanet(nplanet::Int64,jmax::Int64,ntrans::Vector{Int64},params::Vector{T}) where T<:Real
  # Need at least two planets!
  @assert(nplanet>=2)
  # The ntrans vectors should have length nplanet:
  @assert(length(ntrans)==nplanet)
  # Define type of ttv array:
  ttv_el_type = eltype(params) == Float64 ? Float64 : Number
  # Need to create an array to store TTVs with maximum length equal to maximum number of transit times of any planet:
  ntransmax = maximum(ntrans)
  ttv = zeros(T,nplanet,ntransmax) 
  # Each planet requires 5 elements in params (mass_ratio,period,trans0,ecosw,esinw):
  @assert(length(params)==5*nplanet)
  @assert(jmax>=1)  
  for iplanet=1:nplanet
  # Each planet should have at least 2 transits:
    @assert(ntrans[iplanet]>=2)
  end
  for iplanet=1:nplanet-1
  # The periods of the planets should be ordered from least to greatest:
    if (params[(iplanet-1)*5+2] >= params[iplanet*5+2])
      return ttv
    end
  end
  # Set up planets planar-planet types for all of the planets:
  #planet = Array{Planet_plane_hk}(nplanet)
  #planet = Array{Any}(nplanet)
  # Loop over pairs of planets to compute pairwise TTVs
  # Loop over inner planets:
  #println("Looping over planets in ttv_nplanet:")
  for iplanet=1:nplanet-1
    # Create a Planet_plane_hk type for the inner planet:
    p1=Planet_plane_hk(params[(iplanet-1)*5+1],params[(iplanet-1)*5+2],params[(iplanet-1)*5+3],params[(iplanet-1)*5+4],params[(iplanet-1)*5+5])
    # Create an array of times for the inner planet:
    n1 = ntrans[iplanet]
    time1 = collect(p1.trans0 .+ range(0,stop=n1-1,length=n1) .* p1.period)
    # Loop over outer planets:
    for jplanet=iplanet+1:nplanet
      # Create a Planet_plane_hk type for the outer planet:
      p2=Planet_plane_hk(params[(jplanet-1)*5+1],params[(jplanet-1)*5+2],params[(jplanet-1)*5+3],params[(jplanet-1)*5+4],params[(jplanet-1)*5+5])
      # Create an array of times for the outer planet:
      n2 = ntrans[jplanet]
      time2 = collect(p2.trans0 .+ range(0,stop=n2-1,length=n2) .* p2.period)
      # Define arrays to hold the TTVs:
      ttv1=zeros(T,n1)
      ttv2=zeros(T,n2)
      # Call the compute_ttv code which implements equation (33) from Agol & Deck (2016):
      #    println("Calling compute_ttv")
      compute_ttv!(jmax,p1,p2,time1,time2,ttv1,ttv2)
      #    println("Finished compute_ttv")
      for i=1:n1
        ttv[iplanet,i] += ttv1[i]
      end
      for i=1:n2
        ttv[jplanet,i] += ttv2[i]
      end
    end
  end
  return ttv
end
