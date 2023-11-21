from sympy import symbols, Poly

x = symbols('x')
f = (35/128)*x**9 - (180/128)*x**7 + (378/128)*x**5 - (420/128)*x**3 + (315/128)*x
g = (46623/1024)*x**9 - (113492/1024)*x**7 + (97015/1024)*x**5 - (34974/1024)*x**3 + (5850/1024)*x

h = f.subs(x, f.subs(x, g.subs(x, g)))

# h1 = h.subs(x, (x + 4)/16)
# h2 = h.subs(x, (x + 1.85)/16)
# h3 = h.subs(x, (x - 3)/16)

# print(P.all_coeffs()[::-1])

# coeffs = [(h1).coeff(x, i) for i in reversed(range(10))]

# print(coeffs)
# # print((h1+h2+h3).as_coefficients_dict(x))
print(h.subs(x, 0))