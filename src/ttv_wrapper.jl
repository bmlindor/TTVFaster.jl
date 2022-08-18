# Adds mean linear ephemeris and pairwise TTVs from TTVFaster to yield transit times.

include("ttv_nplanet.jl")

function ttv_wrapper(tt0::Vector{Float64},nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer) where T<:Real
  # These lines need modification for different choices of parameters:
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
  for i=1:n2
    tt2[i] += ttv[2,i]
  end
  # If transit times of additional planets were observable these would need to be added in.
  return [tt1;tt2]  
end


#   # If transit times of additional planets were observable these would need to be added in.
#   Currently doesnt account for skipped transits?
# # Arguments:
# - `nplanet::Int64`: number of planets to model
# - `ntrans::Vector{Int64}`: number of transits of each planet in model
# - `params::Vector{Real}`: vector of values for 5 elements of each planet: mass_ratio,period,trans0,ecosw,esinw 
# - `jmax::Int64`: 
# - `treat_earth_moon_as_planet::Bool`: whether to treat the Earth-Moon barycenter as 
# a planet (such that the second planet had no moon)
# # Returns:
# -The model transit times for the observed planets with model TTVs for the system.
function ttv_wrapper(tt0::Vector{Float64},nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer,treat_earth_moon_as_planet::Bool) where T<:Real
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
  for i=1:n2
    if treat_earth_moon_as_planet
      tt2[i] += ttv[2,i]
    else
      # Compute the amplitude of Moon perturbing Earth's transit times:
      tsinphi0 = params[end-2] #tmax sinphi0
      tcosphi0 = params[end-1] #tmax cosphi0
      deltaphi = params[end]
      tt2[i] += ttv[2,i] + tsinphi0*cos((i-1)*deltaphi) + tcosphi0*sin((i-1)*deltaphi)
    end
  end
  return [tt1;tt2]  
end
# function chisquare(tt0,nplanet,ntrans,params,tt,sigtt,jmax,treat_earth_moon_as_planet)
#   chisq = 0.0  #check mearth_moon_as_planetory allocation >>>>>>>>>>>>
#   # println(params,tt[1],sigtt[1])
#   tt_model = ttv_wrapper(tt0,nplanet,ntrans,params,jmax,treat_earth_moon_as_planet) #,fixp3,p3_cur)
#   for j=1:length(tt)
#     chisq += (tt[j]-tt_model[j])^2/sigtt[j]^2
#   end
#   # println(nplanet)
#   return chisq
# end