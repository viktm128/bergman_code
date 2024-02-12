# Want to build program to create g polynomials from N polynomials

using Polynomials
using Plots
include("root_finder.jl")

function build_binomial_expansion(k)
	# Create triangular array from row 1 to m inclusive (omit row 0)
	pt = [[1]]
	for j in 1:k
		new_row = ones(2j + 1)
		for l in 2:2j
			if l % 2 == 0
				new_row[l] = 0
			else
				new_row[l] = last(pt)[l - 2] + last(pt)[l]
			end
		end
		push!(pt, new_row)
	end
	
	return pt
end


function g_coefficients(m,n)
	N = coeffs(rational_coeffs(m,n))
	pt = build_binomial_expansion(m - n)
	new_coeffs = []
	while length(N) > 0
		push!(new_coeffs, N[1])
		N -= N[1]pop!(pt)
		N = N[2:length(N) - 1]
	end

	return new_coeffs
end
