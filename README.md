Why does interpolation work for neural networks? I hope to answer this question by looking at polynomial linear regression with the polynomial degree always being larger than the data.

--What do numerical tests tell us about finite sample behaviour of such functions?

Numerical tests tell us that polynomial interpolation using the least squares solution is numerically terrible. In practice, for high N,P, the polynomial does not really interpolate.

The MP solution gives minimum norm of the solution vector, but this vector is just one of many representations of the polynomial, which may be changed by a change of basis. Different choices of basis give different results. Are any of them linked to other "size" measures, such as the Lp norm of the derivative?


--What are the asymptotics of polynomial interpolation as degree grows faster than the data?

Clearly, we don't have uniform convergence of functions to the truth assuming any amount of noise. Is there convergence is some other sense? Perhaps in l2 of the series expansion of the target function to that of the interpolant?
