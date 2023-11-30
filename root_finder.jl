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

function rational_coeffs(m, n)
	if (gcd(m,n) == 1)
	#	for k in 0:2m - 2n
	#		print("k: ")
	#		print(k)
	#		print(" --> ")
	#		for b2 in 0:2n
	#			print("(")
	#			print(D_1(m,n,k + 2n - 1, b2))
	#			print(",")
	#			print(D_2(m,n,k + 2n - 1, b2))
	#			print(")")
	#		end
	#		println("")
	#	end
		poly = Polynomial([sum([D_1(m,n,k + 2n - 1,b2) * D_2(m,n,k + 2n - 1,b2) for b2 in 0:2n]) for k in 0:2m - 2n], :s)
	end

end

