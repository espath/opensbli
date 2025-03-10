""" Author: djl 05/2022. """
from sympy import flatten, Idx, sqrt, Rational, pprint, factor, nsimplify, collect
from opensbli.core.opensbliobjects import ConstantObject, ConstantIndexed, Globalvariable
from opensbli.core.grid import GridVariable
from opensbli.equation_types.opensbliequations import OpenSBLIEq
from opensbli.core.kernel import Kernel
from opensbli.core.datatypes import Int
from opensbli.core.parsing import EinsteinEquation


class Divergence(object):
    def __init__(self, conservative, inviscid):
        """ Standard divergence form of the equations. Not suitable for most applications, included for educational purposes."""
        self.conservative = conservative
        self.inviscid = inviscid
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        return

    def continuity(self):
        mass = "Eq(Der(rho, t), - Conservative(%s_j, x_j))" % self.rhou
        return mass

    def momentum(self):
        convective_momentum = "(Conservative(rhou_j*u_i, x_j) + Der(p, x_i))"
        return convective_momentum

    def energy(self):
        if self.conservative:
            convective = "(Conservative(rhoE*u_j, x_j))"
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j))" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy

class Blaisdell(object):
    def __init__(self, conservative, inviscid):
        """ Blaisdell spit form - not kinetic energy preserving. Implemented from: 
        G.A. Blaisdell, E.T. Spyropoulos, J.H. Qin. The effect of the formulation of nonlinear terms on aliasing errors in spectral methods. Applied Numerical Mathematics Volume 21, Issue 3, July 1996, Pages 207-219."""
        self.conservative = conservative
        self.inviscid = inviscid
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        return

    def continuity(self):
        mass = "Eq(Der(rho, t), - Skew(rho*u_j, x_j))"
        return mass

    def momentum(self):
        convective_momentum = "Skew(rhou_i*u_j, x_j)"
        return convective_momentum

    def energy(self):
        if self.conservative:
            convective = "(Skew(rhoE*u_j,x_j) - Conservative(p*u_j,x_j))"
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s)" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy


class Jameson(object):
    def __init__(self, conservative, inviscid):
        """ Jameson split form, implemented from:
        A. Jameson. Formulation of Kinetic Energy Preserving Conservative Schemes for Gas Dynamics and Direct Numerical Simulation of One-Dimensional Viscous Compressible Flow in a Shock Tube Using Entropy and Kinetic Energy Preserving Schemes. Journal of Scientific Computing. Volume 34, pages 188–208, (2008)."""
        self.conservative = conservative
        self.inviscid = inviscid
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        return

    def continuity(self):
        mass = "Eq(Der(rho, t), - Conservative(%s_j, x_j))" % self.rhou
        return mass

    def momentum(self):
        convective_momentum = "(1/2) * (Conservative(rhou_i*u_j, x_j) + u_i*Conservative(rhou_j,x_j) + rhou_j * Der(u_i,x_j))"
        return convective_momentum

    def energy(self):
        if self.conservative:
            convective = "(1/2) * (Conservative(rhou_j*H, x_j) + H*Conservative(rhou_j, x_j) + rhou_j*Conservative(H, x_j))"
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j))" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy


class Feiereisen(object):
    def __init__(self, conservative, inviscid):
        """ Quadratic kinetic energy preserving plit form introduced by:
        W. J. Feiereisen, Numerical simulation of a compressible, homogeneous, turbulent shear flow, PhD thesis, Stanford University, 1981."""
        self.conservative = conservative
        self.inviscid = inviscid
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        return

    def continuity(self):
        mass = "Eq(Der(rho, t), - Conservative(%s_j, x_j))" % self.rhou
        return mass

    def momentum(self):
        convective_momentum = "(1/2) * (Conservative(%s_i*u_j, x_j) + %s_j*Der(u_i,x_j) + u_i * Der(%s_j,x_j))" % (self.rhou, self.rhou, self.rhou)
        return convective_momentum

    def energy(self):
        if self.conservative:
            convective = "(1/2) * (Conservative(%s*u_j, x_j) + %s_j*Conservative(%s / rho, x_j) + (%s / rho) * Conservative(%s_j, x_j))" % (self.energy_lhs, self.rhou, self.energy_lhs, self.energy_lhs, self.rhou)
        else:
            convective = "(1/2) * (Conservative(rho*%s*u_j, x_j) + %s_j*Conservative(%s, x_j) + %s * Conservative(%s_j, x_j))" % (self.energy_lhs, self.rhou, self.energy_lhs, self.energy_lhs, self.rhou)
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j))" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s - Conservative(p*u_j, x_j) + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy

class Kok(object):
    def __init__(self, conservative, inviscid):
        """ Kok split form, implemented from:"""
        self.conservative = conservative
        self.inviscid = inviscid
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        return

    def continuity(self):
        mass = "Eq(Der(rho, t), - Conservative(%s_j, x_j))" % self.rhou
        return mass

    def momentum(self):
        convective_momentum = "(1/2) * (Conservative(rhou_i*u_j, x_j) + u_i*Conservative(rhou_j,x_j) + rhou_j * Der(u_i,x_j))"
        return convective_momentum

    def energy(self):
        if self.conservative:
            convective = "( (u_i/2)*Conservative(rhou_i*u_j, x_j) + (rhou_i*u_j/2)*Der(u_i, x_j) + (1/2)*(Conservative(rhou_j*e, x_j) + e*Conservative(rhou_j, x_j) + rhou_j*Der(e, x_j)) + (u_j*Der(p, x_j) + p*Der(u_j, x_j)))"
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s)" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy

class KGP(object):
    def __init__(self, conservative, energy_formulation, inviscid):
        """ KGP split form, implemented with enthalpy splitting from:
        G. Coppola, F. Capuano, S. Pirozzoli, L. de Luca. Numerically stable formulations of convective terms for turbulent compressible flows. Journal of Computational Physics Volume 382, 1 April 2019, Pages 86-104."""
        self.conservative = conservative
        self.inviscid = inviscid
        self.energy_formulation = energy_formulation
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        # KGP coefficients
        self.alpha = Rational(1,4)
        self.beta = Rational(1,4)
        self.delta = Rational(1,4)
        self.gamma = Rational(1,4)
        self.epsilon = 0
        return

    def continuity(self):
        A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
        if self.conservative:
            mass = "Eq(Der(rho, t), - (%s*Conservative(rhou_j, x_j) + %s*Conservative(rhou_j, x_j) + %s*(u_j*Der(rho, x_j) + rho*Der(u_j, x_j)) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
        else:
            mass = "Eq(Der(rho, t), - (%s*Conservative(rho*u_j, x_j) + %s*Conservative(rho*u_j, x_j) + %s*(u_j*Der(rho, x_j) + rho*Der(u_j, x_j)) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
        return mass

    def momentum(self):
        A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
        if self.conservative:
            convective_momentum = "%s*Conservative(rhou_j*u_i, x_j) + %s*(u_i*Conservative(rhou_j, x_j) + rhou_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rhou_i, x_j) + rhou_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
        else:
            convective_momentum = "%s*Conservative(rho*u_j*u_i, x_j) + %s*(u_i*Conservative(rho*u_j, x_j) + rho*u_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rho*u_i, x_j) + rho*u_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
        return convective_momentum

    def energy(self):
        # Split on phi = E, with quadratic split applied to pressure-velocity term
        A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
        if self.energy_formulation == 'enthalpy': # The RHS does not have rhoE or E explicitly here. Pressure divergence derivative is included within H definition H = E + p / rho (constituent relations)
            if self.conservative:
                convective = "(%s*Conservative(rhou_j*H, x_j) + %s*(H*Conservative(rhou_j, x_j) + rhou_j*Conservative(H, x_j)) + %s*(u_j*Conservative(rho*H, x_j) + rho*H*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*H, x_j) + u_j*H*Der(rho, x_j)))" % (A, B, C, D) 
            else:
                convective = "(%s*Conservative(rho*H*u_j, x_j) + %s*(H*Conservative(rho*u_j, x_j) + rho*u_j*Conservative(H, x_j)) + %s*(u_j*Conservative(rho*H, x_j) + rho*H*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*H, x_j) + u_j*H*Der(rho, x_j)))" % (A, B, C, D)
        else:
            if self.conservative:
                convective = "((1/2)*(Conservative(p*u_j, x_j) + p*Der(u_j, x_j) + u_j*Der(p, x_j)) + %s*Conservative(rhoE*u_j, x_j) + %s*((rhoE/rho)*Conservative(rhou_j, x_j) + rhou_j*Conservative((rhoE/rho), x_j)) + %s*(u_j*Conservative(rhoE, x_j) + rhoE*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*(rhoE/rho), x_j) + u_j*(rhoE/rho)*Der(rho, x_j)))" % (A, B, C, D)
            else:
                convective = "((1/2)*(Conservative(p*u_j, x_j) + p*Der(u_j, x_j) + u_j*Der(p, x_j)) + %s*Conservative(rho*Et*u_j, x_j) + %s*(Et*Conservative(rho*u_j, x_j) + rho*u_j*Conservative(Et, x_j)) + %s*(u_j*Conservative(rho*Et, x_j) + rho*Et*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*Et, x_j) + u_j*Et*Der(rho, x_j)))" % (A, B, C, D)
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s)" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy


class KEEP(object):
    def __init__(self, conservative, energy_formulation, inviscid, pressure_fix=True):
        """ KEEP scheme implemented from:
        Y Kuya, K Totani, S Kawai. Kinetic energy and entropy preserving schemes for compressible flows by split convective forms. Journal of Computational Physics, (2018). """
        self.conservative = conservative
        self.inviscid = inviscid
        self.pressure_fix = True
        self.energy_formulation = energy_formulation
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        # KGP coefficients
        self.alpha = Rational(1,4)
        self.beta = Rational(1,4)
        self.delta = Rational(1,4)
        self.gamma = Rational(1,4)
        self.epsilon = 0
        return

    def continuity(self):
        A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
        if self.conservative:
            mass = "Eq(Der(rho, t), - (%s*Conservative(rhou_j, x_j) + %s*Conservative(rhou_j, x_j) + %s*(u_j*Der(rho, x_j) + rho*Der(u_j, x_j)) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
        else:
            mass = "Eq(Der(rho, t), - (%s*Conservative(rho*u_j, x_j) + %s*Conservative(rho*u_j, x_j) + %s*(u_j*Der(rho, x_j) + rho*Der(u_j, x_j)) + %s*(rho*Der(u_j, x_j) + u_j*Der(rho, x_j))))" % (A, B, C, D)
        return mass

    def momentum(self):
        A, B, C, D = self.alpha, self.beta, self.gamma, self.delta
        if self.conservative:
            convective_momentum = "%s*Conservative(rhou_j*u_i, x_j) + %s*(u_i*Conservative(rhou_j, x_j) + rhou_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rhou_i, x_j) + rhou_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
        else:
            convective_momentum = "%s*Conservative(rho*u_j*u_i, x_j) + %s*(u_i*Conservative(rho*u_j, x_j) + rho*u_j*Der(u_i, x_j)) + %s*(u_j*Conservative(rho*u_i, x_j) + rho*u_i*Der(u_j, x_j)) + %s*(rho*Conservative(u_j*u_i, x_j) + u_i*u_j*Der(rho, x_j))" % (A, B, C, D)
        return convective_momentum

    def energy(self):
        if self.pressure_fix:
            A, B = Rational(1,2), Rational(1,2)
            convective = "(%s*((u_k/2)*Conservative(rhou_k*u_j, x_j) + (rhou_k*u_j/2)*Der(u_k, x_j) + (rhou_k/2)*Conservative(u_k*u_j, x_j) + (u_k*u_j/2)*Conservative(rhou_k, x_j)) \
            + %s * (Conservative(rhou_j*e, x_j) + u_j*Conservative(rho*e, x_j) + rho*e*Der(u_j, x_j)) \
            + u_j*Der(p, x_j) + p*Der(u_j, x_j))" % (A, B)
        else:
            A, B = Rational(1,2), Rational(1,4)
            convective = "(%s*((u_k/2)*Conservative(rhou_k*u_j, x_j) + (rhou_k*u_j/2)*Der(u_k, x_j) + (rhou_k/2)*Conservative(u_k*u_j, x_j) + (u_k*u_j/2)*Conservative(rhou_k, x_j)) \
            + %s * (Conservative(rhou_j*e, x_j) + u_j*Conservative(rho*e, x_j) + e*Conservative(rhou_j, x_j) + rho*Conservative(e*u_j, x_j)+ rho*e*Der(u_j, x_j) + rho*u_j*Der(e, x_j) + e*u_j*Der(rho, x_j)) \
            + u_j*Der(p, x_j) + p*Der(u_j, x_j))" % (A, B)
        if self.inviscid:
            energy = "Eq(Der(%s, t), - %s)" % (self.energy_lhs, convective)
        else:
            energy = "Eq(Der(%s, t), - %s + Der(q_j, x_j) + Der(u_i*tau_i_j, x_j))" % (self.energy_lhs, convective)
        return energy


class NS_Split(object):
    """ Split forms for the convective parts of the Navier-Stokes equations with central/DRP schemes."""
    def __init__(self, split_type, ndim, constants, coordinate_symbol="x", conservative=True, viscosity=None, energy_formulation='none', debug=False):
        self.viscosity = viscosity
        # Add diffusive terms?
        if self.viscosity == 'inviscid':
            self.inviscid = True
        else:
            self.inviscid = False
        # Which splitting method to use
        if split_type == 'Divergence':
            print("Convective terms are using the straight Divergence form.")
            self.split = Divergence(conservative, self.inviscid)
        elif split_type == 'Blaisdell':
            print("Convective terms are using the Blaisdell split form.")
            self.split = Blaisdell(conservative, self.inviscid)
        elif split_type == 'Feiereisen':
            print("Convective terms are using the Feiereisen split form.")
            self.split = Feiereisen(conservative, self.inviscid)
        elif split_type == 'Jameson':
            print("Convective terms are using the Jameson split form.")
            self.split = Jameson(conservative, self.inviscid)
        elif split_type == 'Kok':
            print("Convective terms are using the Kok split form.")
            self.split = Kok(conservative, self.inviscid)
        elif split_type == 'KGP':
            print("Convective terms are using the Kennedy-Gruber-Pirozzoli split form.")
            self.split = KGP(conservative, energy_formulation, self.inviscid)
        elif split_type == 'KEEP':
            print("Convective terms are using the KEEP split form.")
            self.split = KEEP(conservative, energy_formulation, self.inviscid)
        else:
            raise NotImplementedError("Only Divergence, Blasidell, Feiereisen, Jameson, Kok, KGP, and KEEP splitting methods are implemented.")

        self.conservative = conservative
        self.coordinate_symbol = coordinate_symbol
        self.constants = constants
        self.ndim = ndim

        self.EE = EinsteinEquation()
        self.replace_factors = False
        # Storing either the conservative or primitive variables as the q vector to advance in time.
        if self.conservative:
            self.rhou = 'rhou'
            self.mom_lhs = 'rhou'
            self.energy_lhs = 'rhoE'
        else:
            self.rhou = 'rho*u'
            self.mom_lhs = 'u'
            self.energy_lhs = 'Et'
        # Viscous and heat-flux substitutions
        if debug: # Don't expand the diffusive terms
            self.substitutions = []
        else:
            if self.inviscid:
                self.substitutions = []
            else:
                self.substitutions = self.diffusive_terms()
        self.mass = self.continuity_eq()
        self.momentum = self.momentum_eq()
        self.energy = self.energy_eq()
        return

    def continuity_eq(self):
        out = self.split.continuity()
        out = self.EE.expand(out, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        return out

    def momentum_eq(self):
        if self.inviscid:
            base_momentum = "Eq(Der(%s_i, t), - Der(p, x_i))" % self.mom_lhs
        else:
            base_momentum = "Eq(Der(%s_i, t), - Der(p, x_i) + Der(tau_i_j, x_j))" % self.mom_lhs
        out = self.EE.expand(base_momentum, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        # Add convective parts
        convective = self.split.momentum()
        expanded_convective = self.EE.expand(convective, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        expanded_convective[0] = factor(expanded_convective[0])
        for no, value in enumerate(out):
            temp = OpenSBLIEq(out[no].lhs,  out[no].rhs - expanded_convective[no])
            out[no] = temp
        return out

    def energy_eq(self):
        energy = self.split.energy()
        out = self.EE.expand(energy, self.ndim, self.coordinate_symbol, self.substitutions, self.constants)
        return out

    def diffusive_terms(self):
        """ Viscous stress tensor and heat-flux terms, depending on whether viscosity is variable or not."""
        if self.viscosity == 'constant':
            stress_tensor = "Eq(tau_i_j, (1.0/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)*Der(u_k,x_k)))"
            heat_flux = "Eq(q_j, ((1.0/Re)/((gama-1)*Minf*Minf*Pr))*Der(T,x_j))"
        else:
            stress_tensor = "Eq(tau_i_j, (mu/Re)*(Der(u_i,x_j)+ Der(u_j,x_i)- (2/3)* KD(_i,_j)*Der(u_k,x_k)))"
            heat_flux = "Eq(q_j, ((mu/Re)/((gama-1)*Minf*Minf*Pr))*Der(T,x_j))"
        substitutions = [stress_tensor, heat_flux]
        return substitutions
