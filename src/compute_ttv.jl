include("ttv_succinct.jl")
"""
Computes transit timing variations to linear order in eccentricity for non-resonant, plane-parallel planets.
Please cite Agol & Deck (2015) if you make use of this in published research.
"""

struct Planet_plane # Deprecated, only exported to make the error message work
  mass_ratio :: Float64
  period   :: Float64
  trans0   :: Float64
  ecc    :: Float64  
  omega    :: Float64   
end

struct Planet_plane_hk{T<:Real} # Parameters of a planet in a plane-parallel system
  # Mass ratio of the planet to the star:
  mass_ratio :: T
  # Initial time of transit:
  period   :: T
  trans0   :: T
  # e times cos or sin of longitude of periastron measured from line of sight, in radians:
  ecosw    :: T
  esinw    :: T
end

"""
compute_ttv(jmax,p1,p2,time1,time2,ttv1,ttv2)

Arguments:
  jmax:  Maximum j over which to sum the TTV calculation for both planets
  p1:    Planet type for inner planet
  p2:    Planet type for outer planet
  time1: Transit times for inner planet
  time2: Transit times for outer planet
  ttv1:  Container for TTVs of the inner planet
  ttv2:  Container for TTVs of the outer planet
Outputs: 
  ttv1
  ttv2
"""

function compute_ttv!(jmax::Integer,p1::Planet_plane_hk{T},p2::Planet_plane_hk{T},time1::Vector{T},time2::Vector{T},ttv1::Vector{T},ttv2::Vector{T}) where T<:Real

  # Compute the semi-major axis ratio of the planets:
  # println(p1.period,p2.period)
  alpha = (p1.period/p2.period)^(2//3)  
  #println(alpha, p1.period, p2.period)
  @assert(alpha < 1)
  @assert(alpha > 0)
  # Number of times:
   ntime1 = length(time1)
   ntime2 = length(time2)
  f1=zeros(T,jmax+2,5)  
  f2=zeros(T,jmax+2,5)  
  # Compute the coefficients:
  ttv_succinct!(jmax+1,alpha,f1,f2)  # I need to compute coefficients one higher than jmax
  # Compute TTVs for inner planet (equation 33):
  # Compute since of \pomegas:
  e1 = sqrt(p1.esinw*p1.esinw + p1.ecosw*p1.ecosw)
  e2 = sqrt(p2.esinw*p2.esinw + p2.ecosw*p2.ecosw)
  sin1om=p1.esinw/e1
  sin2om=p2.esinw/e2
  cos1om=p1.ecosw/e1
  cos2om=p2.ecosw/e2
  # Compute mean motions:
  n1=2pi/p1.period
  n2=2pi/p2.period
  # Compute initial longitudes:
  lam10=-n1*p1.trans0 + 2*p1.esinw # 2*p1.eccen*sin1om
  lam20=-n2*p2.trans0 + 2*p2.esinw # 2*p2.eccen*sin2om
  @inbounds for i=1:ntime1
  # Compute the longitudes of the planets at times of transit of planet 1 (equation 49):
    lam11 = n1*time1[i]+lam10
    lam21 = n2*time1[i]+lam20
    psi1  = lam11-lam21 # Compute difference in longitudes at times of transit of planet 1
    sinpsi1=sin(psi1)
    cospsi1=cos(psi1)
    sinlam11 = sin(lam11)
    coslam11 = cos(lam11)
    sinlam1om1=sinlam11*cos1om - coslam11*sin1om
    coslam1om1=coslam11*cos1om + sinlam11*sin1om
    sinlam1om2=sinlam11*cos2om - coslam11*sin2om
    coslam1om2=coslam11*cos2om + sinlam11*sin2om
    ttv1[i]=zero(T) #0.0
    sinjm1psi1=zero(T) #0.0
    cosjm1psi1=one(T) #1.0
  # Sum over j:
    for j=1:jmax
      sinjpsi1=sinjm1psi1*cospsi1 + cosjm1psi1*sinpsi1
      cosjpsi1=cosjm1psi1*cospsi1 - sinjm1psi1*sinpsi1
      ttv1[i] += f1[j+1,1]*sinjpsi1
      ttv1[i] += f1[j+1,2]*e1*(sinjpsi1*coslam1om1 - cosjpsi1*sinlam1om1)
      ttv1[i] += f1[j+1,3]*e1*(sinjpsi1*coslam1om1 + cosjpsi1*sinlam1om1)
      ttv1[i] += f1[j  ,4]*e2*(sinjpsi1*coslam1om2 - cosjpsi1*sinlam1om2)
      ttv1[i] += f1[j+2,5]*e2*(sinjpsi1*coslam1om2 + cosjpsi1*sinlam1om2)
      sinjm1psi1=sinjpsi1
      cosjm1psi1=cosjpsi1
    end
  # Multiply by period and mass ratio, and divide by 2*Pi:
    ttv1[i] = ttv1[i]*p1.period*p2.mass_ratio/(2pi)
  end
  # Compute TTVs for outer planet (equation 33):
  @inbounds for i=1:ntime2
  # Compute the longitudes of the planets at times of transit of planet 2:
    lam12 = n1*time2[i]+lam10
    lam22 = n2*time2[i]+lam20
    sinlam22 = sin(lam22)
    coslam22 = cos(lam22)
    psi2  = lam12-lam22 # Compute difference in longitudes at times of transit of planet 2
    sinpsi2=sin(psi2)
    cospsi2=cos(psi2)
    sinlam2om1=sinlam22*cos1om - coslam22*sin1om
    coslam2om1=coslam22*cos1om + sinlam22*sin1om
    sinlam2om2=sinlam22*cos2om - coslam22*sin2om
    coslam2om2=coslam22*cos2om + sinlam22*sin2om
    ttv2[i]=zero(p2.period) #0.0
    sinjm1psi2=zero(p2.period) #0.0
    cosjm1psi2=one(p2.period) #1.0
  # Sum over j:
    for j=1:jmax
      sinjpsi2=sinjm1psi2*cospsi2 + cosjm1psi2*sinpsi2
      cosjpsi2=cosjm1psi2*cospsi2 - sinjm1psi2*sinpsi2
      ttv2[i] += f2[j+1,1]*sinjpsi2
      ttv2[i] += f2[j+1,2]*e2*(sinjpsi2*coslam2om2 - cosjpsi2*sinlam2om2)
      ttv2[i] += f2[j+1,3]*e2*(sinjpsi2*coslam2om2 + cosjpsi2*sinlam2om2)
      ttv2[i] += f2[j+2,4]*e1*(sinjpsi2*coslam2om1 - cosjpsi2*sinlam2om1)
      ttv2[i] += f2[j  ,5]*e1*(sinjpsi2*coslam2om1 + cosjpsi2*sinlam2om1)
      sinjm1psi2=sinjpsi2
      cosjm1psi2=cosjpsi2
    end
  # Multiply by period and mass ratio, and divide by 2*Pi:
    ttv2[i] = ttv2[i]*p2.period*p1.mass_ratio/(2pi)
  end
  # Finished! 
  return 
end  # compute_ttv!

function compute_ttv!(jmax::Integer,p1::Planet_plane,p2::Planet_plane,time1::Vector,time2::Vector,ttv1::Vector,ttv2::Vector)
  # Error message to explain to anyone who tries to use the old version
  error("The Planet_plane data structure has been deprecated in favor of Planet_plane_hk")
end