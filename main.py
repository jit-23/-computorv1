import sys
import re as re
import cmath
from collections import defaultdict

PI = 3.14159265358979323846

def my_cos(x):
    while x > PI: x -= 2 * PI
    while x < -PI: x += 2 * PI
    term = 1.0
    total = 1.0
    for n in range(1, 20):
        term *= -x * x / ((2 * n - 1) * (2 * n))
        total += term
    return total

def my_acos(x):
    lo, hi = 0.0, PI
    for _ in range(100):
        mid = (lo + hi) / 2
        if my_cos(mid) > x:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def my_sqrt(x):
    last_guess= x/2.0
    while True:
        guess= (last_guess + x/last_guess)/2
        if abs(guess - last_guess) < .000001: # example threshold
            return guess
        last_guess= guess


def cbrt(x):
    return -((-x) ** (1/3)) if x < 0 else x ** (1/3)

def main():
	if len(sys.argv)  < 2:
		Polinomial = input("")
	else:
		Polinomial = sys.argv[1]
	if Polinomial == "":
		sys.stderr.write("Error: Invalid Polymonial format\n")
		return
	if Polinomial.count('=') != 1:
		sys.stderr.write("Error: equation must contain exactly one '='\n")
		return
	left_eq, right_eq = Polinomial.lower().replace(" ","").split('=')
	
	term_re = re.compile(r'([+-]?\d+(?:\.\d+)?)\*x\^(\d+)')
	full_term_re = re.compile(r'^(?:[+-]?\d+(?:\.\d+)?\*x\^\d+)+$')
	
	def is_valid_side(s):
		return s == "0" or full_term_re.fullmatch(s)

	if not is_valid_side(left_eq) or not is_valid_side(right_eq):
		sys.stderr.write("Error: Invalid Polymonial format\n")
		exit(1)

	left_terms = re.findall(term_re, left_eq)
	right_terms = re.findall(term_re, right_eq)
	

	terms = [(float(c), int(p)) for c, p in left_terms]
	terms += [(-float(c), int(p)) for c, p in right_terms]

	coeffs = defaultdict(float)      # it behaves just  like std::map<int, double>
	for c, p in terms:
		coeffs[p] += c     # -> it will put the terms in order of power to x and the '+c' will sum all the bases on the same power


	#---

	print("Reduced form:", end=" ")

	eq = []

	for i , p in enumerate(sorted(coeffs, reverse=True)):
		c = coeffs[p]
		sign = "+" if c >= 0 else "-"
		if i == 0:
			eq.append(f"{c if c >= 0 else -c} * X^{p}" if c >= 0 else f"-{-c} * X^{p}")
		else:
			eq.append(f"{sign} {abs(c)} * X^{p}")
	print(" ".join(eq), "= 0")



	degree = max(coeffs, default=0)
	while degree > 0 and coeffs[degree] == 0:
		del coeffs[degree]           # a dict has no "last element", so remove by key
		degree -= 1

	print(f"Polinomial degree: {degree}")
	#	quadratic formula: ax^2 + bx + c = 0
	if degree == 3:
		a,b,c,d = coeffs[3], coeffs[2], coeffs[1], coeffs[0]
		p = (3*a*c - b**2)/(3*a**2)
		q = (27 * a**2 * d + 2 * b**3 - 9*a*b*c)/(27 * a**3)

		delta = -(4 * p**3) - (27 * q**2)
		# de solving of hela diferent depending of the delta signal(-/+)
		if delta < 0:
			D = (q/2)**2 + (p/3)**3
			u = cbrt(-(q/2) + (my_sqrt(D)))
			v = cbrt(-(q/2) - (my_sqrt(D)))
			x1 = u + v - b / (3*a) ## x1 = y - b/3*a || y = u + v
			x2 = -(u + v)/2 - b / (3*a) + (u - v)* my_sqrt(3)/(2j)
			x3 = -(u + v)/2 - b / (3*a) - (u - v)* my_sqrt(3)/(2j)
			print("There is 1 real root and 2 complex roots:")
			print(f"x1 = {x1.real:.6f}")
			print(f"x2 = {x2:.6f}")
			print(f"x3 = {x3:.6f}")
		elif delta == 0:
			u = cbrt((-q/2))
			v = u 
			x1 =2 * u - b / (3*a)
			x2 = x3 = -u - b / (3*a)
			print("There are 3 real roots, two are equal:")
			print(f"x1 = {x1.real:.6f}")
			print(f"x2 = {x2.real:.6f}")
			print(f"x3 = {x3.real:.6f}")
		else :
			r = my_sqrt(-p / 3)
			arg = (3 * q) / (2 * p) * my_sqrt(-3 / p)
			arg = max(-1.0, min(1.0, arg))
			theta = my_acos(arg)
			roots = [2 * r * my_cos((theta - 2 * PI * k) / 3) - b / (3*a) for k in range(3)]
			x1, x2, x3 = roots
			print("There are 3 distinct real roots:")
			print(f"x1 = {x1:.6f}")
			print(f"x2 = {x2:.6f}")
			print(f"x3 = {x3:.6f}")

			#D = (q/2)**2 + (p/3)**3
			#u = (-(q/2) + my_sqrt(D))**(1/3)
			#v = (-(q/2) - my_sqrt(D))**(1/3)
			#x1 = u + v - (b / (3*a))
			#x2 = -((u + v)/2 - b / (3*a) + (u - v))* my_sqrt(3)/(2j)
			#x3 = -((u + v)/2 - b / (3*a) - (u - v))* my_sqrt(3)/(2j)
			#print("There are 3 distinct real roots:")
			#print(f"x1 = {x1.real:.6f}")
			#print(f"x2 = {x2.real:.6f}")
			#print(f"x3 = {x3.real:.6f}")

	elif degree == 2:
		a, b, c = coeffs[2], coeffs[1], coeffs[0]      # by power: order and missing terms don't matter
		delta = b ** 2 - 4 * a * c

		if (delta > 0):
			sqr = my_sqrt(delta)
			solution1 = (-b + sqr) / (2 * a)
			solution1+=0.0
			solution2 = (-b - sqr) / (2 * a)
			solution2+=0.0
			
			print("Discriminant is strictly positive, the two solutions are:");print(f"{solution1:.6f}");print(f"{solution2:.6f}")
		elif (delta == 0):
			solution = -b / (2 * a)
			solution+=0.0
			print("Discriminant is zero, the solution is:");print(f"{solution:.6f}")
		else:
			print("Discriminant is strictly negative, the two complex solutions are:")
			sqr = my_sqrt(-delta)
			real_part = -b / (2 * a)
			real_part+=0.0
			imaginary_part = sqr / (2 * a)
			imaginary_part+=0.0
			print(f"{real_part:.6f} + {imaginary_part:.6f}i")
			print(f"{real_part:.6f} - {imaginary_part:.6f}i")
	elif degree == 1:
		b,c = coeffs[1], coeffs[0]
		solution = -c / b
		solution+=0.0
		print("The solution is:");print(f"{solution:.6f}")
	
	elif degree == 0:
		if coeffs[0] == 0:
			print("All real numbers are solutions.")
		else:
			print("No solution exists.")
	else:
		print("The polynomial degree is greater than 2, I can't solve.")

main()
