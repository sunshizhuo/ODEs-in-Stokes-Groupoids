import sympy
from sympy import exp, Inverse, Rational, log

class Sto_k():
    # The kth Stokes Groupoid
    # Over the space of C acting on C
    def __init__(self, k):
        assert(type(k) == int)
        assert(k >= 1)
        self.k = k
    def source(self, u, z):
        return z
    def target(self, u, z):
        return (exp(u * (z**(self.k-1))) * z).simplify()
    def product(self, p2, p1):
        # Compute p1 \circ p2
        # p for path.
        u1, z1 = p1
        u2, z2 = p2
        if self.source(u2, z2) == self.target(u1, z1):
            ansu = (u2*exp((self.k-1)*u1*z1**(self.k-1)) + u1).simplify()
            return (ansu, z1)
        else:
            raise("The product is not defined")
    def inverse(self, u, z):
        ansu = (-u * exp(-(self.k-1)*u*z**(self.k-1))).simplify()
        return (ansu, self.target(u, z))

def get_phi_from_psi(psi, k):
    from sympy.abc import u, z
    S = Sto_k(k)
    psi_end = psi.subs(z, S.target(u, z))
    psi_start = psi.subs(z, S.source(u, z))
    prod = psi_end * Inverse(psi_start)
    prod = prod.simplify()
    return prod

def solve_constant_A(k, A):
    from sympy.abc import z
    # assert(A.subs(z, 0) == A)
    # In this case, A is constant.
    if k == 0:
        return exp(A*z)
    if k == 1:
        return exp(A*log(z))
    else:
        return exp(Rational(1, 1-k) * A * z**(1-k))
        
def get_Phi_from_simple_ode(k, A):
    psi = solve_constant_A(k, A)
    return get_phi_from_psi(psi, k)