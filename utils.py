import sympy
from sympy import Derivative, Function, Symbol, symbols, Eq, pi, cos, sin, exp, log, oo
from sympy import Function, dsolve, Derivative, checkodesol
from sympy import fps, Rational
from sympy import pprint, Matrix, eye, zeros
from sympy import Inverse
from MatrixSeries import *
from Stokes import *

psi = Function('psi')
phi = Function('phi')
z = Symbol('z')

def ignore_const(u, val):
    for i in u.free_symbols: 
        if i != z:
            u = u.subs(i, val)
    return u

def solve_1d(k, poly, psi):
    eq = Eq(z**k * Derivative(psi(z), z), poly * psi(z))
    result = dsolve(eq, psi(z))
    return result.simplify()

def solve_diag(k, diag_entries):
    psi = Function('psi')
    psis = []
    n = len(diag_entries)
    for i in range(0, n):
        poly = diag_entries[i]
        psii = solve_1d(k, poly, psi)
        ## We do this since we already know that psi is diagonal.
        psii = ignore_const(psii.rhs, 1)
        psis.append(psii)
    return Matrix(psis)

def solve_phi(k, A):
    AA = MatrixSeriesA(A)
    n = AA.n
    B = sympy.zeros(1, n)
    for i in range(0, k):
        B += AA[i].diagonal() * z**i
    psi = solve_diag(k, B)
    Psi = sympy.diag(list(psi), unpack=True)
    Phi = get_phi_from_psi(Psi, k)
    return Psi, Phi


def get_order(A, off_diag):
    """
    order = get_order(A, False)
    off_diag_order = get_order(A, True)
    # Hp are all 0 after off_diag_order terms
    """
    order = 0
    n = A.shape[0]
    for i in range(0, n):
        for j in range(0, n):
            if off_diag and (i == j): continue
            if A[i, j].is_polynomial():
                order = sympy.Max(sympy.degree(A[i, j], z), order)
            else:
                order = oo
    return order


def get_Hp(AA, k, p):
    """
    if Aps[i, j] == 0 or i == j:
    Note that if any entry is infinite after the division, this means that such entry could be any number. 
    This is because we just need to satisfy (alpha_i - alpha_j)Hp[i, j] = Aps[i, j]
    so it is infinite iff alpha_i = alpha_j. In this case, Hp[i, j] could be any number.
    """
    assert(k >= 0)
    Aps = AA[p]
    n = AA.n
    if k == 1:
        return Matrix(n, n, lambda i, j: 0 if i == j else Aps[i, j]/(AA.alphas[i]-AA.alphas[j]-p))
    if k > 1:
        return Matrix(n, n, lambda i, j: 0 if i == j else Aps[i, j]/(AA.alphas[i]-AA.alphas[j]))

def transform_up_to(Series, order):
    l = order
    Ans = zeros(Series.n)
    for i in range(0, l):
        Ans += Series[i]*z**i
    return Ans

def Gauge(k, F, A):
    F_inverse = F.inverse()
    F_B = (F * A + Derivative(F, z).simplify() * z**k) * F_inverse 
    F_B = F_B.simplify()
    return F_B

def get_diagonal(A):
    A_diag = A.diagonal()
    ans = sympy.diag(list(A_diag), unpack=True)
    return ans


def get_Gauge_up_to_order(A, k, order):
    AA = MatrixSeriesA(A)
    f = lambda i: get_Hp(AA, k, i)
    F = construct_from_Hp(AA.n, z, f)
    ## This is the first Gauge Transform
    
    A_diag = get_diagonal(A)
    if A_diag.free_symbols == set():
        ## it is constant
        K_trunc = eye(AA.n)
    else:
        AA_diag = MatrixSeriesA(A_diag)
        F_after = truncate_after_k(AA_diag, k)
        ## SS_first_k = truncate_first_k(SS, k)
        log_K = integrate(F_after)
        Ans = transform_up_to(log_K, order)
        K_trunc = exp(-Ans).simplify()
        ## This is the truncated second Gauge Transform
    F_trunc = transform_up_to(F, order)
    Total = K_trunc * F_trunc
    Total = sympy.simplify(Total)
    ## The final Gauge transform is KS, corresponding to the transform of K \circ S
    return K_trunc, F_trunc, Total




