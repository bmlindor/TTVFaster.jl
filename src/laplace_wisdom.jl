"""
 # Compute Laplace coefficients and Leverrier derivative/
 #         j
 #    j   d     i
 #   a   ---   b (a)
 #         j    s
 #       da

 #  by series summation 
 # Code due to Jack Wisdom 
"""
function laplace_wisdom(s::Rational,i::Integer,j::Integer,a::T) where T<:Real

 # define LAPLACE_EPS 1.0e-12
  LAPLACE_EPS = convert(T,1.0e-12)

  #if (i lt 0) then i = -i
  i=abs(i)

  if j <= i     # compute first term in sum 
    factor4 = one(T)
    for k=0:j-1
      factor4 *= i - k
    end
    lap_coef_sum = factor4
    q0 = 0
  else
    q0 = fld(j + 1 - i,2) # largest integer less than or equal to x/y
    lap_coef_sum = zero(T)
    factor4 = one(T)
  end

  # Compute factors for terms in lap_coef_sum: 
  factor1 = s
  factor2 = s + i
  factor3 = i + 1
  for q=1:q0-1       # no contribution for q = 0 
    factor1 *= s + q
    factor2 *= s + i + q
    factor3 *= i + q + 1
  end
  if q0 > 1
    q=q0
  else
    q=1
  end
  #println(j+1-i,q0)
  term = a*a * factor1 * factor2 / (factor3 * q)

  # Sum series: 
  while  (term*factor4) > LAPLACE_EPS
    factor4 = one(T)
    for k=0:j-1
      factor4 *= (2*q + i - k)
    end
    lap_coef_sum += term * factor4
    factor1 += 1
    factor2 += 1
    factor3 += 1
    q = q+1
    term *= a*a * factor1 * factor2 / (factor3 * q)
  end

  # Fix coefficient:
  for k=0:i-1
    lap_coef_sum *= (s+k)/(k+1)
  end

  apower = (q0 <= 0) ?  i : 2*q0 + i - 2
  lap_coef_sum *= 2 * a^apower
  # Return the Laplace Coefficient:
  return lap_coef_sum
end