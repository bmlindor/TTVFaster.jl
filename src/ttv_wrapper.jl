# Adds pairwise TTVs to mean linear ephemeris to calculate transit times?
# this doesnt account for skipped transits

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
  ttv1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) 
  for i=1:n1
    ttv1[i]+= ttv[1,i]
  end
  t02 = params[8]
  per2 = params[7]
  ttv2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2)) 
  for i=1:n2
    ttv2[i] += ttv[2,i]
  end
  # If transit times of additional planets were observable these would need to be added in.
  #println("param2: ",param)
  return [ttv1;ttv2]  
end

function ttv_wrapper(tt0::Vector{Float64},nplanet::Int64,ntrans::Vector{Int64},params::Vector{T},jmax::Integer,treat_earth_moon_as_planet::Bool) where T<:Real
  # These lines need modification for different choices of parameters:
  n1,n2 = ntrans[1:2]
  # println(n1, " ", n2)

  # Call ttv_nplanet:
  ttv = ttv_nplanet(nplanet,jmax,ntrans,params[1:5*nplanet])
  # We measure transit times,not TTVs,so add back in the linear ephemeris:
  t01 = params[3]
  per1 = params[2]
  ttv1 = collect(range(t01,stop = t01+per1*(n1-1),length = n1)) 
  for i=1:n1
    ttv1[i]+= ttv[1,i]
  end
  t02 = params[8]
  per2 = params[7]
  ttv2 = collect(range(t02,stop = t02+per2*(n2-1),length = n2)) 
  for i=1:n2
    if treat_earth_moon_as_planet
      ttv2[i] += ttv[2,i]
    else
      # Compute the amplitude of Moon perturbing Earth's transit times:
      tsinphi0 = params[end-2] #tmax sinphi0
      tcosphi0 = params[end-1] #tmax cosphi0
      deltaphi = params[end]
      ttv2[i] += ttv[2,i] + tsinphi0*cos((i-1)*deltaphi) + tcosphi0*sin((i-1)*deltaphi)
    end
  end
  # If transit times of additional planets were observable these would need to be added in.
  #println("param2: ",param)
  return [ttv1;ttv2]  
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