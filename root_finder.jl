using Polynomials
using Plots

function generate_integer_poly(k)
	s = Polynomial([0, 1], :s)
	p_k = Polynomial([n * (k - n) for n in 1:(k - 1)], :s)
	q_k = Polynomial([n^2 for n in 1:k], :s)
	q_k += Polynomial([(k - n)^2 for n in 1:k], :s) * s^k 
	# do not need to multiply by an extra s from the t term because we will factor it out
	p_k * s + q_k + p_k * s^(k - 1)  # reduce the multiplication factor on the p terms for the same reason
end

theta = range(0, 2*pi, length=100)
x = cos.(theta)
y = sin.(theta)

plot(x, y, lw=2, lc=:black)


for k = 2:10
	r = roots(generate_integer_poly(k))
	scatter!(real(r), imag(r), label="k = "*string(k)) 
end

savefig("bergman_roots_2-10.png")
