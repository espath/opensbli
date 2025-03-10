"""@brief
   @author David J Lusher
   @contributors Satya Pramod Jammy
   @details
"""
from sympy import diag, eye, Rational, pprint, simplify
from opensbli.core.opensbliobjects import ConstantObject, EinsteinTerm, DataSet
from sympy.parsing.sympy_parser import parse_expr


class EulerEquations(object):
    """ Class to generate the Eigensystems used to diagonalize the Euler equations."""

    def __init__(self, ndim, species=None, **kwargs):
        self.ndim = ndim
        self.species = species
        return

    def apply_direction(self, direction):
        """ Substitutes a direction to get the eigensystem for that direction.

        :arg int direction: Direction to obtain the eigensystem."""
        # Dictionaries to store ev, REV and LEV for each direction
        ev_dict, LEV_dict, REV_dict = {}, {}, {}
        # Metric symbols from block
        met_symbols = self.met_symbols
        n_metric_terms = len(met_symbols.atoms(DataSet))
        # Metric terms for this direction to substitute into the matrix
        terms = [EinsteinTerm('k%d' % i) for i in range(self.ndim)]
        if n_metric_terms > 0:
            metric_values = [self.detJ*met_symbols[direction, i] for i in range(self.ndim)]
        else:
            metric_values = [met_symbols[direction, i] for i in range(self.ndim)]
        subs_dict = dict([(x, y) for (x, y) in zip(terms, metric_values)])
        # Scaling factor based on metrics
        if n_metric_terms > 0:
            factor = self.detJ*sum([met_symbols[direction, i]**2 for i in range(self.ndim)])**(Rational(1, 2))
        else:
            factor = sum([met_symbols[direction, i]**2 for i in range(self.ndim)])**(Rational(1, 2))
        required_metrics = factor.atoms(DataSet)
        subs_dict[EinsteinTerm('k')] = factor

        def g(x):
            return x.subs(subs_dict, evaluate=False)
        ev_dict[direction] = diag(*list(self.ev.applyfunc(g)))
        LEV_dict[direction] = self.LEV.applyfunc(g)
        REV_dict[direction] = self.REV.applyfunc(g)
        if n_metric_terms > 0:
            # remove the detJ from the 1 / sqrt(D00*2 + D10**2) factor, as it cancels out in the LEV/REV matrices
            factor = factor / self.detJ
        return ev_dict, LEV_dict, REV_dict, required_metrics, factor

    def generate_eig_system(self, block):
        """ Creates the Eigensystems used to diagonalize the Euler equations. No direction is
        passed at this stage, the Eigensystems are in general form.

        :arg object block: OpenSBLI SimulationBlock.
        :returns: None """
        ndim = self.ndim
        # Check if block has metrics
        metrics = block.fd_metrics

        if metrics == eye(block.ndim):
            self.met_symbols = eye(block.ndim)
            self.detJ = 1
        elif metrics.is_diagonal() and metrics != eye(block.ndim):
            self.met_symbols = eye(block.ndim)
            self.detJ = block.detJ_metrics[0]
        else:
            self.met_symbols = metrics
            self.detJ = block.detJ_metrics[0]
        local_dict = {'Symbol': EinsteinTerm, 'gama': ConstantObject('gama')}

        if ndim == 1:
            matrix_symbols = ['H']
            matrix_formulae = ['a**2/(gama-1) + u0**2/2']
            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            if self.species == 'passive_scalar':
                REV = 'Matrix([[-2.0/a**2, 0, 2.0/(a**2*(gama - 1.0)), 2.0/(a**2*(gama - 1.0))], [-2.0*u0/a**2, 0, 2.0*(-a + u0)/(a**2*(gama - 1)), 2.0*(a + u0)/(a**2*(gama - 1))], [-1.0*u0**2/a**2, 0, 1.0*(2.0*a**2 - 2.0*a*gama*u0 + 2.0*a*u0 + 1.0*gama*u0**2 - 1.0*u0**2)/(a**2*(1.0*gama**2 - 2.0*gama + 1.0)), 1.0*(2.0*a**2 + 2.0*a*gama*u0 - 2.0*a*u0 + 1.0*gama*u0**2 - 1.0*u0**2)/(a**2*(1.0*gama**2 - 2.0*gama + 1.0))], [-2.0*f/a**2, 1.0/a, 2.0*f/(a**2*(gama - 1.0)), 2.0*f/(a**2*(gama - 1.0))]])'
                LEV = 'Matrix([[-0.5*a**2 + 0.25*gama*u0**2 - 0.25*u0**2, 0.5*u0*(-gama + 1), 0.5*gama - 0.5, 0], [-1.0*a*f, 0, 0, 1.0*a], [0.25*u0*(a + 0.5*u0*(gama - 1.0))*(gama - 1.0), -0.25*(a + u0*(gama - 1.0))*(gama - 1.0), 0.25*(gama - 1.0)**2, 0], [0.25*u0*(-a + 0.5*u0*(gama - 1.0))*(gama - 1.0), 0.25*(a - u0*(gama - 1.0))*(gama - 1.0), 0.25*(gama - 1.0)**2, 0]])'
                ev = 'diag([u0, u0, -a + u0, a + u0])'
            elif self.species == 'N N2':
                REV = 'Matrix([[-1, -2*(rhoN + rhoN2)/a, 2*rhoN/(a*(gama - 1)), 2*rhoN/(a*(gama - 1))], [1, 0, 2*rhoN2/(a*(gama - 1)), 2*rhoN2/(a*(gama - 1))], [0, -2*u0*(rhoN + rhoN2)/a, 2*(-a + u0)*(gama - 1)*(rhoN + rhoN2)/(a*(gama - 1)**2), 2*(a + u0)*(rhoN + rhoN2)/(a*(gama - 1))], [0, -1*u0**2*(rhoN + rhoN2)/a, (rhoN + rhoN2)*(2*a**2*(gama - 1) + 2*a**2 - 2*a*gama*u0*(gama - 1.0) + 1.0*gama*u0**2*(gama - 1.0))/(a*gama*(gama - 1.0)**2), (rhoN + rhoN2)*(2.0*a**2*(gama - 1.0) + 2.0*a**2 + 2.0*a*gama*u0*(gama - 1.0) + 1.0*gama*u0**2*(gama - 1.0))/(a*gama*(gama - 1.0)**2)]])'
                LEV = 'Matrix([[-0.5*rhoN2*u0**2*(gama - 1)/(a**2*(rhoN + rhoN2)), (1.0*a**2*rhoN + 1.0*a**2*rhoN2 - 0.5*gama*rhoN2*u0**2 + 0.5*rhoN2*u0**2)/(a**2*(rhoN + rhoN2)), 1.0*rhoN2*u0*(gama - 1)/(a**2*(rhoN + rhoN2)), -1.0*rhoN2*(gama - 1)/(a**2*(rhoN + rhoN2))], [1.0*(-0.5*a**2 + 0.25*gama*u0**2 - 0.25*u0**2)/(a*(rhoN + rhoN2)), 1.0*(-0.5*a**2 + 0.25*gama*u0**2 - 0.25*u0**2)/(a*(rhoN + rhoN2)), -0.5*u0*(gama - 1)/(a*(rhoN + rhoN2)), 0.5*(gama - 1)/(a*(rhoN + rhoN2))], [0.25*u0*(a + 0.5*u0*(gama - 1.0))*(gama - 1.0)/(a*(rhoN + rhoN2)), 0.25*u0*(a + 0.5*u0*(gama - 1.0))*(gama - 1.0)/(a*(rhoN + rhoN2)), 0.25*(a*(-gama + 1) - u0*(gama - 1.0)**2)/(a*(rhoN + rhoN2)), 0.25*(gama - 1.0)**2/(a*(rhoN + rhoN2))], [0.25*u0*(-a + 0.5*u0*(gama - 1.0))*(gama - 1.0)/(a*(rhoN + rhoN2)), 0.25*u0*(-a + 0.5*u0*(gama - 1.0))*(gama - 1.0)/(a*(rhoN + rhoN2)), 0.25*(a*(gama - 1) - u0*(gama - 1.0)**2)/(a*(rhoN + rhoN2)), 0.25*(gama - 1.0)**2/(a*(rhoN + rhoN2))]])'
                ev = 'diag([u0, u0, -a + u0, a + u0])'
            else:
                ev = 'diag([u0-a, u0, u0+a])'
                REV = 'Matrix([[1,1,1], [u0-a,u0,u0+a], [H-u0*a,u0**2 /2,H+u0*a]])'
                LEV = 'Matrix([[ u0*(2*H + a*u0 - u0**2)/(2*a*(2*H - u0**2)), (-H - a*u0 + u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)],[2*(H - u0**2)/(2*H - u0**2),2*u0/(2*H - u0**2), 2/(-2*H + u0**2)],[u0*(-2*H + a*u0 + u0**2)/(2*a*(2*H - u0**2)),  (H - a*u0 - u0**2/2)/(a*(2*H - u0**2)),  1/(2*H - u0**2)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)
            pprint(REV)
            print('------------------------------------------------')
            pprint(LEV)
            print('------------------------------------------------')

            subs_dict = dict(zip(matrix_symbols, matrix_formulae))

            def f(x):
                return x.subs(subs_dict)
            self.ev = ev.applyfunc(f)
            self.REV = REV.applyfunc(f)
            self.LEV = LEV.applyfunc(f)

        if ndim == 2:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2.0))', '1/(rho*a*sqrt(2.0))', 'k0*u0 + k1*u1', '(gama-1)*((u0**2 + u1**2)/2)', 'sqrt(k0**2 + k1**2)']
            matrix_symbols_2 = ['k0_t', 'k1_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k0_t*u0 + k1_t*u1', 'k0*u0+k1*u1']

            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[1,0,alpha,alpha], [u0,k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [u1,-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [phi_sq/(gama-1),rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + a*theta_t),alpha*((phi_sq+a**2)/(gama-1) - a*theta_t)]])'
            LEV = 'Matrix([[1-phi_sq/a**2,(gama-1)*u0/a**2,(gama-1)*u1/a**2,-(gama-1)/a**2], [-(k1_t*u0-k0_t*u1)/rho,k1_t/rho,-k0_t/rho,0], [bta*(phi_sq-a*theta_t),bta*(k0_t*a-(gama-1)*u0),bta*(k1_t*a-(gama-1)*u1),bta*(gama-1)], [bta*(phi_sq+a*theta_t),-bta*(k0_t*a+(gama-1)*u0),-bta*(k1_t*a+(gama-1)*u1),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)

            # Apply the sub
            def f(x):
                return x.subs(subs_dict, evaluate=False)

            def g(x):
                return x.subs(subs_dict2, evaluate=False)
            self.ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            self.REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            self.LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)

        elif ndim == 3:
            matrix_symbols = ['alpha', 'bta', 'theta', 'phi_sq', 'k']
            matrix_formulae = ['rho/(a*sqrt(2.0))', '1/(rho*a*sqrt(2.0))', 'k0*u0 + k1*u1 + k2*u2', '(gama-1)*((u0**2 + u1**2 + u2**2)/2)', 'sqrt(k0**2 + k1**2 + k2**2)']
            matrix_symbols_2 = ['k0_t', 'k1_t', 'k2_t', 'theta_t', 'U']
            matrix_formulae_2 = ['k0/k', 'k1/k', 'k2/k', 'k0_t*u0 + k1_t*u1 + k2_t*u2', 'k0*u0+k1*u1+k2*u2']

            matrix_symbols = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols]
            matrix_formulae = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae]
            matrix_symbols_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_symbols_2]
            matrix_formulae_2 = [parse_expr(l, local_dict=local_dict, evaluate=False) for l in matrix_formulae_2]
            subs_dict = dict(zip(matrix_symbols, matrix_formulae))
            subs_dict2 = dict(zip(matrix_symbols_2, matrix_formulae_2))

            ev = 'diag([U, U, U, U+a*k, U-a*k])'
            REV = 'Matrix([[k0_t,k1_t,k2_t,alpha,alpha], [k0_t*u0,k1_t*u0-k2_t*rho,k2_t*u0+k1_t*rho,alpha*(u0+k0_t*a),alpha*(u0-k0_t*a)], [k0_t*u1+k2_t*rho,k1_t*u1,k2_t*u1-k0_t*rho,alpha*(u1+k1_t*a),alpha*(u1-k1_t*a)], [k0_t*u2-k1_t*rho,k1_t*u2+k0_t*rho,k2_t*u2,alpha*(u2+k2_t*a),alpha*(u2-k2_t*a)], [k0_t*phi_sq/(gama-1) + rho*(k2_t*u1-k1_t*u2),k1_t*phi_sq/(gama-1) + rho*(k0_t*u2-k2_t*u0),k2_t*phi_sq/(gama-1) + rho*(k1_t*u0-k0_t*u1),alpha*((phi_sq+a**2)/(gama-1) + theta_t*a),alpha*((phi_sq+a**2)/(gama-1) - theta_t*a)]])'
            LEV = 'Matrix([[k0_t*(1-phi_sq/a**2)-(k2_t*u1-k1_t*u2)/rho,k0_t*(gama-1)*u0/a**2,k0_t*(gama-1)*u1/a**2 + k2_t/rho,k0_t*(gama-1)*u2/a**2 - k1_t/rho,-k0_t*(gama-1)/a**2], [k1_t*(1-phi_sq/a**2)-(k0_t*u2-k2_t*u0)/rho,k1_t*(gama-1)*u0/a**2 - k2_t/rho,k1_t*(gama-1)*u1/a**2,k1_t*(gama-1)*u2/a**2 + k0_t/rho,-k1_t*(gama-1)/a**2], [k2_t*(1-phi_sq/a**2) - (k1_t*u0-k0_t*u1)/rho,k2_t*(gama-1)*u0/a**2 + k1_t/rho,k2_t*(gama-1)*u1/a**2 - k0_t/rho,k2_t*(gama-1)*u2/a**2,-k2_t*(gama-1)/a**2], [bta*(phi_sq-theta_t*a),-bta*((gama-1)*u0-k0_t*a),-bta*((gama-1)*u1-k1_t*a),-bta*((gama-1)*u2 - k2_t*a),bta*(gama-1)], [bta*(phi_sq+theta_t*a),-bta*((gama-1)*u0+k0_t*a),-bta*((gama-1)*u1+k1_t*a),-bta*((gama-1)*u2+k2_t*a),bta*(gama-1)]])'
            ev = parse_expr(ev, local_dict=local_dict, evaluate=False)
            REV = parse_expr(REV, local_dict=local_dict, evaluate=False)
            LEV = parse_expr(LEV, local_dict=local_dict, evaluate=False)

            # Apply the sub
            def f(x):
                return x.subs(subs_dict, evaluate=False)

            def g(x):
                return x.subs(subs_dict2, evaluate=False)
            self.ev = ev.applyfunc(f).applyfunc(g).applyfunc(g)
            self.REV = REV.applyfunc(f).applyfunc(g).applyfunc(g)
            self.LEV = LEV.applyfunc(f).applyfunc(g).applyfunc(g)
        return
