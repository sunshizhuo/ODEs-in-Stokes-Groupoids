import sympy
from sympy import Derivative, Function, Symbol, symbols, Eq, pi, cos, sin, exp, log, oo
from sympy import Function, dsolve, Derivative, simplify
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

def derivative(F, z, l):
    """
    Misteriously, Derivative sometimes return constant 0.
    """
    n = F.shape[0]
    ans = Derivative(F, z, l)
    if F == 0:
        ans = zeros(n)
    return ans

def solve_phi(k, G):
    if G.free_symbols == set():
        B = G
    else:
        B = G_up_to(G, k-1)
        # print(G)
        # AA = MatrixSeriesA(G)
        # n = AA.n
        # B = sympy.zeros(1, n)
        # for i in range(0, k):
        #     B += AA[i].diagonal() * z**i
    B = B.diagonal()
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


def get_F(AA, k, p):
    """
    Find the F from 1 to p (inclusive) inductively.
    Note that if any entry is infinite after the division, this means that such entry could be any number. 
    This is because we just need to satisfy the generic condition.
    """
    assert(k >= 1)
    n = AA.n
    A = AA.value
    F = eye(n)
    for l in range(1, p+1):
        if(A.free_symbols == set()):
            break
        Ais = (Derivative(A, z, l).subs(z, 0) / sympy.factorial(l)).simplify()
        assert(Ais.is_Matrix and Ais.shape == (n,n))
        if k == 1:
            His = Matrix(n, n, lambda i, j: 0 if i == j else Ais[i, j]/(AA.alphas[i]-AA.alphas[j]-l))
        else:
            His = Matrix(n, n, lambda i, j: 0 if i == j else Ais[i, j]/(AA.alphas[i]-AA.alphas[j]))
        Fi = eye(n) + His * z**l
        F = F + His*F*(z**l)
        A = Gauge(k, Fi, A)
    return F, A

def transform_up_to(Series, order):
    l = order
    Ans = zeros(Series.n)
    for i in range(0, l):
        Ans += Series[i]*z**i
    return Ans

def Gauge(k, F, A):
    F_inverse = F.inv()
    F_B = (F * A + simplify(F.diff(z)) * z**k) * F_inverse
    F_B = simplify(F_B)
    return F_B

def get_g(k, A):
    n = A.shape[0]
    A0 = simplify(A.subs(z, 0))
    A0_diag = A0.diagonal()
    g = [0] * n
    eigenvals = set()
    for i in range(0, n):
        while A0_diag[i]+g[i] in eigenvals:
            g[i]+=1
        eigenvals.add(A0_diag[i]+g[i])
        g[i] = solve_constant_A(k, g[i])

    final_g = sympy.diag(g, unpack=True)
    return final_g

def get_diagonal(A):
    A_diag = A.diagonal()
    ans = sympy.diag(list(A_diag), unpack=True)
    return ans

def G_up_to(G, ord):
    M = G
    G_approx = simplify(G.subs(z, 0))
    for i in range(1, ord+1):
        M = Derivative(M, z) / i
        coeff = simplify(M.subs(z,0))
        G_approx += coeff * (z**i)
    return G_approx



def get_Gauge_up_to_order(A, k, order):
    AA = MatrixSeriesA(A)
    # f = lambda i: get_Hp(AA, k, i)
    # F = construct_from_Hp(AA.n, z, f)
    F_trunc, G = get_F(AA, k, order)
    ## This is the first Gauge Transform
    if G.free_symbols == set():
        K_trunc = eye(AA.n)
    else:
        G_approx = G_up_to(G, order)
        G_approx = MatrixSeriesA(G_approx)
        F_after = truncate_after_k(G_approx, k)
        ## SS_first_k = truncate_first_k(SS, k)
        log_K = integrate(F_after)
        Ans = transform_up_to(log_K, order)
        K_trunc = exp(-Ans).simplify()
        ## This is the truncated second Gauge Transform
    # F_trunc = transform_up_to(F, order)
    Total = K_trunc * F_trunc
    Total = sympy.simplify(Total)
    ## The final Gauge transform is KS, corresponding to the transform of K \circ S
    return K_trunc, F_trunc, Total




