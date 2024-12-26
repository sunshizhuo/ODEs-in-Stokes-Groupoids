import sympy
from sympy import Matrix, fps, Rational


def show(expr):
    from IPython.display import Latex, display, Math
    latex_str = "$"+sympy.latex(expr)+"$"
    display(Math(latex_str))


class MatrixSeries():
    def __init__(self, n, z):
        self.n = n
        self.variable = z
    def get_diagonalD(self, i):
        Aps = self[i]
        return Aps.diagonal()

class MatrixSeriesFromGenerator(MatrixSeries):
    def __init__(self, n, z, generator, start):
        super().__init__(n, z)
        self.generator = generator
        self.start = start

    def __getitem__(self, k):
        if self.start <= k:
            return self.generator(k)
        else:
            return sympy.zeros(self.n)

    def __str__(self):
        z = self.variable
        i = self.start
        expr = self[i]*z**i + self[i+1]*z**(i+1) + self[i+2]*z**(i+2)
        show(expr)
        return "$"+sympy.latex(expr)+"...........$"
    

class MatrixSeriesA(MatrixSeries):
    def __init__(self, A):
        assert(A.shape[0] == A.shape[1])
        super().__init__(A.shape[0], A.free_symbols.pop())
        self.value = A
        
        A0 = A.subs(self.variable, 0)
        assert(A0.is_diagonal())
        self.alphas = A0.diagonal()
        
        self.series = Matrix(self.n, self.n, lambda i, j: fps(A[i, j]))
        

    def __getitem__(self, k):
        z = self.variable
        def entry(i, j):
            # try:
            #     return self.series[i, j][k] / (z**k)
            # except:
            #     ## In this case, Aseries[i, j] is a constant
            #     return self.series[i, j] if k == 0 else 0
            if self.value[i, j].is_constant():
                ## In this case, Aseries[i, j] is a constant
                return self.value[i, j] if k == 0 else 0
            else:
                return self.series[i, j][k] / (z**k)
        return Matrix(self.n, self.n, entry)


def integrate(Series):
    n, z = Series.n, Series.variable
    return MatrixSeriesFromGenerator(n, z, (lambda i: (Rational(1, i) * Series[i-1]) if i >= 1 else zeros(n)), Series.start+1)

def truncate_first_k(Series, k):
    n, z = Series.n, Series.variable
    Ans = sympy.zeros(n)
    for i in range(Series.start, k):
        Ans += Series[i]*z**i
    return Ans

def truncate_after_k(Series, k):
    n, z = Series.n, Series.variable
    return MatrixSeriesFromGenerator(n, z, (lambda i: Series[i+k]), 0)


def construct_from_Hp(n, z, get_nth_factor):
    ## This is the series from the product prod_{i=1}^infty(I+z^iH_i)
    ## Each coefficient is H_k + H_1H_{k-1} + H_2H_{k-2} + ...
    def generate(k):
        if k < 0:
            return sympy.zeros(n)
        if k == 0:
            return sympy.eye(n)
        ans = get_nth_factor(k)
        for i in range(1, (k+1)//2):
            ans += get_nth_factor(i) * get_nth_factor(k-i)
        return ans
    NewSeries = MatrixSeriesFromGenerator(n, z, generate, 0)
    return NewSeries
        