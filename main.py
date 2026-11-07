from sys import argv
import re as re
from collections import defaultdict

def my_sqrt(x):
    if x == 0:
        return 0.0
    g = x if x > 1 else 1.0
    for _ in range(100):
        g = (g + x / g) / 2
    return g


def main():
	Polinomial = argv[1]		
	left_eq, right_eq = Polinomial.lower().replace(" ","").split('=')
	
	term_re = re.compile(r'([+-]?\d+(?:\.\d+)?)\*x\^(\d+)')

	left_terms = re.findall(term_re, left_eq)
	right_terms = re.findall(term_re, right_eq)
	

	terms = [(float(c), int(p)) for c, p in left_terms]
	terms += [(-float(c), int(p)) for c, p in right_terms]

	print(terms)
	print("=")
	coeffs = defaultdict(float)      # it behaves just  like std::map<int, double>
	for c, p in terms:
		coeffs[p] += c     # -> it will put the terms in order of power to x and the '+c' will sum all the bases on the same power
	for p in coeffs:
		print(p, coeffs[p])

	#---

	degree = max(coeffs, default=0)
	while degree > 0 and coeffs[degree] == 0:
		del coeffs[degree]           # a dict has no "last element", so remove by key
		degree -= 1

	#	quadratic formula: ax^2 + bx + c = 0

	if degree == 2:
		a, b, c = coeffs[2], coeffs[1], coeffs[0]      # by power: order and missing terms don't matter
		delta = b ** 2 - 4 * a * c

		if (delta > 0):
			sqr = my_sqrt(delta)
			solution1 = (-b + sqr) / (2 * a)
			solution2 = (-b - sqr) / (2 * a)
			print("Discriminant is strictly positive, the two solutions are:");print(solution1);print(solution2)
		elif (delta == 0):
			solution = -b / (2 * a)
			print("Discriminant is zero, the solution is:");print(solution)
		else:
			print("Discriminant is strictly negative, the two complex solutions are:")
			sqr = my_sqrt(-delta)
			real_part = -b / (2 * a)
			imaginary_part = sqr / (2 * a)
			print(f"{real_part} + {imaginary_part}i")
			print(f"{real_part} - {imaginary_part}i")
	elif degree == 1:
		b,c = coeffs[1], coeffs[0]
		solution = -c / b
		print("The solution is:");print(solution)
	
	elif degree == 0:
		if coeffs[0] == 0:
			print("All real numbers are solutions.")
		else:
			print("No solution exists.")
	else:
		print("The polynomial degree is greater than 2, I can't solve.")

## amanha faco o 3grau sozinho

main()
