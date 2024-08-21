using Polynomials
using Plots

function generate_integer_poly(k, theta=0)
	s = Polynomial([0, 1], :s)
	p_k = Polynomial([n * (k - n) for n in 1:(k - 1)], :s)
	q_k = Polynomial([n^2 for n in 1:k], :s)
	q_k += Polynomial([(k - n)^2 for n in 1:k], :s) * s^k 

	# do not need to multiply by an extra s from the t term because we will factor it out
	p_k * exp(im * theta) * s + q_k * exp(im * theta) + p_k * s^(k - 1)  # reduce the multiplication factor on the p terms for the same reason
end

function up_to_K_max(K_MAX)
	theta = range(0, 2*pi, length=100)
	x = cos.(theta)
	y = sin.(theta)
	
	plot(x, y, lw=2, lc=:black)


	for k = 2:10
		r = roots(generate_integer_poly(k))
		scatter!(real(r), imag(r), label="k = "*string(k)) 
	end

	savefig("roots_2-10.png")

end


function rotate_with_theta(k)
	theta = range(0, (2 - 2 / 100) * pi, 100)
	r = roots.(generate_integer_poly.(k, theta))
	
	psi = range(0, 2*pi, 100)
	x = cos.(psi)
	y = sin.(psi)

	plot(x, y, lw=2, lc=:red, legend = false, size=(1200,800))
	scatter!(real.(r), imag.(r), mc=:black, legend = false)
	#xlims!(minimum(real.(r)) - 1, maximum(real.(r)) + 1)
	#ylims!(minimum(imag.(r)) - 1, maximum(imag.(r)) + 1)
	savefig("k=" * string(k) * "_rotational_roots.png")
end

function D_1(m, n, d, b_2)
	b_1 = d - b_2
	if (0 <= b_1 <= m - 2)
		val = b_1 + 1
	elseif (m - 1 <= b_1 <= 2m - 2)
		val = 2m - b_1 - 1
	else
		val = 0
	end
end

function D_2(m, n, d, b_2)
	alpha = 2m*n - n
	gamma = (m - n)b_2 + n*d
	if (alpha - m + 1 <= gamma <= alpha)
		val = -alpha + m + gamma
	elseif (alpha + 1 <= gamma <= alpha + m - 1)
		val = alpha + m - gamma
	else
		val = 0
	end
end

# technically, there is a factor of n thats been pulled out already (hidden in the D1 value of each term)
function rational_coeffs(m, n)
	@assert m > n "m needs to be greater than n"
	@assert gcd(m, n) == 1 "m/n must be a reduced fraction"

	poly = Polynomial([sum([D_1(m,n,k + 2n - 1,b2) * D_2(m,n,k + 2n - 1,b2) for b2 in 0:2n]) for k in 0:2m - 2n], :s)

end


function d1_array(m,n)
	# finds all coefficients involved in sum of rational_coeffs polynomial evaluated at s=-1
	# summing rowwise across output should match the coeffs of rational_coeffs up to a negative sign
	@assert m > n "m needs to be greater than n"
	@assert gcd(m, n) == 1 "m/n must be a reduced fraction"
	d1 = [D_1(m, n, d, b2) for b2 in 0:2n for d in (2n - 1):(2m-1)]
	d1 = reshape(d1, 2(m - n) + 1, 2n + 1)
end

function d2_array(m,n)
	# finds all coefficients involved in sum of rational_coeffs polynomial evaluated at s=-1
	# summing rowwise across output should match the coeffs of rational_coeffs up to a negative sign
	@assert m > n "m needs to be greater than n"
	@assert gcd(m, n) == 1 "m/n must be a reduced fraction"
	d2 = [D_2(m, n, d, b2) for b2 in 0:2n for d in (2n - 1):(2m-1)]
	d2 = reshape(d2, 2(m - n) + 1, 2n + 1)
end

function summands_array(m, n)
  d1_array(m,n) .* d2_array(m, n)
end



function diagonal_sum(m,n)
	# test intution that evaluating the summands array by adding on the counter diagonals would make it easier
	# failed intution - counter example --> m/n = 7/2
	a = summands_array(m,n)
	b = []
	for j in 0:2m
		push!(b, 0)
		if j <= 2n
			startrow = 0
			startcol = j + 2
		else
			startrow = j - 2n
			startcol = 2n + 2
		end

		for t in 1:min(j + 1, 2m - 2n + 1, 2m - j + 1, 2n + 1)
			b[length(b)] += a[startrow + t, startcol  - t]
		end
	end

	return b
end


function intuition_check_diagonals(M_MAX, N_MAX)
	# iterate diagonal_sum across some collection of m/n to test for counter examples
	for n in 2:N_MAX
		for m in (n+1):M_MAX
			if gcd(m,n) == 1
				if any(diagonal_sum(m,n) .> 0)
				       print([m,n])
				end
			end
		end
	end
end

function monotonicity_intuition_check(d, ITER_MAX)
	@assert d % 2 > 0 "d must be odd"
	
	prev = 0
	for n = 2:ITER_MAX
		m = n + d
		if gcd(m,n) == 1
			new = rational_coeffs(m,n)(-1)
			if new > prev
				print([m,n])
			else
				prev = new
			end
		end
	end
	println("All cases passed")
	
end

function cosine_coeffs(m,n)
  rc = [sum([D_1(m,n,k + 2n - 1,b2) * D_2(m,n,k + 2n - 1,b2) for b2 in 0:2n]) for k in 0:2m - 2n]
  return reverse(rc[1:m-n+1])
end

function fejer_differences(m,n)
  rc = [sum([D_1(m,n,k + 2n - 1,b2) * D_2(m,n,k + 2n - 1,b2) for b2 in 0:2n]) for k in 0:2m - 2n]
  fc = reverse(rc[1:m-n+1])
  return fc - push!(fc[2:m-n+1], 0)
end


function sub_E(m, n, j)
  return floor(Int, ((j + 1) * n - 1) // m)
end

function sub_f(m, n, j)
  return (j + 1) * SparsePolynomial(Dict(n=>1), :s) + (m - j - 1) * SparsePolynomial(Dict(m=>1), :s)
end

function sub_g(m, n, j)
  ej = sub_E(m, n, j)
  return Polynomial(Int.(n * [(j + 1 - m // n * ej), (m // n + m // n * ej - j - 1)]), :s)
end

function sub_kernels(m, n)
  return [sub_f(m, n, j) * sub_g(m, n, j) * SparsePolynomial(Dict((j + n - 1 - sub_E(m, n, j))=>1), :s) for j in 0:m-1]
end
