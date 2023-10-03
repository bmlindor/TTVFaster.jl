include("ttv_nplanet.jl")
""" 
ttv_wrapper(tt0)

Adds mean linear ephemeris and pairwise TTVs from TTVFaster to yield transit times. 
Currently doesnt account for skipped transits?
Arguments:
  tt0:     Initial times of transit for planets
  nplanet: Number of planets to model
  ntrans:  Number of transits for each planet in model
  params:  Holds 5 elements of each planet: mass_ratio,period,trans0,ecosw,esinw 
  jmax:    Maximum j over which to sum the TTV calculation for both planets
Returns:
  tts:  The model transit times for the observed planets with model TTVs for the system.
"""
function ttv_wrapper(tt0::Vector{Float64},nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer) where T<:Real
  # Call ttv_nplanet:
  ttv = ttv_nplanet(nplanet,jmax,ntrans,params[1:5*nplanet])
  # We measure transit times,not TTVs,so add back in the linear ephemeris:
  tts=Float64[] 
    for iplanet in 1:nplanet
      t0=params[(iplanet-1)*5+3]
      per=params[(iplanet-1)*5+2]
      n=ntrans[iplanet]
      tt=collect(range(t0,stop=t0+per*(n-1),length=n))
      # time1=collect(t0.+ range(0,stop=n-1,length=n) .*per)
      # println(time1[1], " VS ", tt[1])
      # println("For planet ",iplanet,", T0=",t0," per=", per," ntrans=",n)
      tt.+=ttv[iplanet,1:n] 
      append!(tts,tt) 
    end
  return tts
end

"""
If nplanet=2 and we want to model that planet 2 hosts a satellite in a circular orbit: 
Arguments (if different from above):
  params: Holds 5 elements of each planet, plus 3 for satellite  
  treat_PMB_as_planet: Whether to treat the Planet-Moon Barycenter as an observed planet for p2 (i.e. for case where Earth has no Moon)

Returns:
  tts:   The model transit times for the observed planets with model TTVs for the system.
"""
function ttv_wrapper(tt0::Vector{Float64},nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer,treat_PMB_as_planet::Bool) where T<:Real
  # If transit times of additional planets were observable these would need to be added in.
  n1,n2 = ntrans[1:2]
  # println(n1, " ", n2)

  # Call ttv_nplanet:
  ttv = ttv_nplanet(nplanet,jmax,ntrans,params[1:5*nplanet])
  # We measure transit times,not TTVs,so add back in the linear ephemeris:
  t01 = params[3]
  per1 = params[2]
  tt1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) 
  for i=1:n1
    tt1[i]+= ttv[1,i]
  end
  t02 = params[8]
  per2 = params[7]
  tt2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2)) 
  # If modeling satellite in circular orbit about another planet, these lines need modification for different choices of parameters.
  for i=1:n2
    if treat_PMB_as_planet
      tt2[i] += ttv[2,i]
    else
    # Compute the TTVs of planet 2, given satellite:
    # Satellite needs 3 elements in params: (tsinphi0,tcosphi0,deltaphi) , where t is maximum amplitude of TTVs on planet 2
      tsinphi0 = params[end-2] 
      tcosphi0 = params[end-1] 
      deltaphi = params[end]
      tt2[i] += ttv[2,i] + tsinphi0*cos((i-1)*deltaphi) + tcosphi0*sin((i-1)*deltaphi)
    end
  end
  return [tt1;tt2]  
end
