from sympy import symbols
from sympy import tan, pi, tanh, sinh, cosh, cos, sin, exp, summation, sqrt
from numpy import random
import numpy as np
from opensbli import *

def forcing_setup(n_terms, direction, block, seed_count):
	evaluations = []
	theta, f, g, h, t = symbols('theta f g h t', **{'cls': GridVariable})
	omega, dt, xa, xb, A = symbols('omega dt xa xb A', **{'cls': ConstantObject})
	# Current iteration
	current_iter = block.get_temporal_schemes[0].iteration_number
	# Coordinates
	x0, x1, x2 = symbols('x0 x1 x2', **{'cls': DataObject})
	Lx0, Lx1, Lx2 = symbols('Lx0 Lx1 Lx2', **{'cls': ConstantObject})
	if direction == 1:
		span_coord, span_len = x2, Lx2
	elif direction == 2:
		span_coord, span_len = x1, Lx1
	# Number of terms
	L, M = n_terms, n_terms
	T1 = GridVariable('T1')
	# Random numbers
	random.seed(2*seed_count)
	l_randoms = [random.uniform(0,1) for _ in range(n_terms)]
	random.seed(2*seed_count+1)
	m_randoms = [random.uniform(0,1) for _ in range(n_terms)]
	# theta_eqn = Eq(theta, 2*pi*(x0 - xa)/(xb - xa))
	# f_eqn = Eq(f, 4*sin(theta)*(1-cos(theta))/sqrt(27.0))
	f_eqn = Eq(f, exp(-(x0-(xb-xa))**2 / 16))

	t_eqn = Eq(T1, (1-Rational(4,5))/(1-Rational(4,5)**M))

	# Corner scaling equation
	# scaling_eqn = corner_scaling(direction, block, span_coord, span_len)
	# evaluations += [scaling_eqn]
	# scale = scaling_eqn.lhs
	time_eqn = Eq(t, dt*current_iter)
	# Manually create g_eqn sum
	g_eqn = 0
	for i in range(1,n_terms+1):
		g_eqn += t_eqn.rhs*Rational(4,5)**(i-1) * sin(2*pi*i*(span_coord/span_len + l_randoms[i-1]))
	g_eqn = Eq(g, g_eqn)

	# Manually create h_eqn sum
	omegas = np.linspace(0.04, 0.12, n_terms)
	np.random.shuffle(omegas)
	# print(omegas)
	h_eqn = 0
	for i in range(1,n_terms+1):
		h_eqn += t_eqn.rhs*Rational(4,5)**(i-1) * sin(omegas[i-1]*t + 2*pi*m_randoms[i-1])
	h_eqn = Eq(h, h_eqn)

	evaluations += [f_eqn, time_eqn, g_eqn, h_eqn]
	# final_eqn = scale*A*f*g*h
	final_eqn = A*f*g*h
	print("Forcing equations")
	for eqn in evaluations:
		pprint(eqn)
	expr_condition_pairs = Piecewise((0.0, x0<xa), (0.0, x0>xb), (final_eqn, True))
	norm_vel_eqn = Eq(GridVariable('u_norm'), expr_condition_pairs)
	evaluations += [norm_vel_eqn]
	return evaluations
