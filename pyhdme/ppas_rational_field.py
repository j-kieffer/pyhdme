r"""
Principally polarized abelian surfaces over QQ and complex period matrices

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

from sage.rings.complex_arb import ComplexBallField
from pyhdme.ppas_invariants import PPASInvariants

# Copyright 2025 Jean Kieffer
# See LICENSE file for license details.

class PPASRationalField(SageObject):
    r"""b

    Create a data structure encoding a principally polarized abelian surface
    over the rational field, allowing computations with period matrices in
    interval arithmetic at varying precisions.
    """

    def __init__(self, data, inv_type = "Modular"):
        r"""
        Initialize a PPASRationalField object. The input data should be as
        specified by the PPASInvariants class, over the field of rational
        numbers.
        """

        prec = 32
        CC = ComplexBallField(prec)

        self.__qq_invariants = PPASInvariants.__init__(data, inv_type)
        self.__cplx_prec = prec
        self.__cplx_invariants = self.__qq_invariants.change_ring(CC, force_exact_computations = True)
        self.__cplx_roots_prec = 0
        self.__cplx_roots = None
        self.__cplx_periods_prec = 0
        self.__cplx_periods = None

    def rational_invariants(self):
        return self.__qq_invariants

    def complex_invariants(self, prec = self.__cplx_prec):
        if self.__cplx_prec < prec:
            self.__cplx_prec = prec
            # this includes base-changing the curve equation if it was already computed.
            self.__cplx_invariants = self.rational_invariants().change_ring(CC)
            return self.__cplx_invariants
        elif self.__cplx_prec > prec:
            CC = ComplexBallField(prec)
            return self.__cplx_invariants.change_ring(CC)
        else:
            return self.__cplx_invariants

    def is_geometrically_split(self):
        return self.rational_invariants().is_geometrically_split()

    def complex_field(self):
        return self.complex_invariants().base_ring()

    def complex_genus_2_curve(self, prec = self.__cplx_prec):
        if self.is_geometrically_split():
            raise ValueError("This PPAS is geometrically split")
        return self.complex_invariants(prec).genus_2_curve_equation()

    def complex_weierstrass_points(self, prec = self.__cplx_prec):
        if self.__cplx_roots is None:
            p = max(prec, 32)
            done = False
            while not done:
                self.__cplx_roots_prec = p
                try:
                    rts = Sequence(self.complex_genus_2_curve(p).roots(multiplicities = False),
                                   universe = self.complex_field())
                except ValueError: #unable to isolate roots
                    #todo: this will happen if the complex curve is a degree 5 model...
                    continue
                if len(rts) < 6:
                    continue
                done = True
                for i in range(7):
                    for j in range(i + 1, 7):
                        if rts[i].overlaps(rts[j]):
                            done = False
                p *= 2
            self.__cplx_roots = rts
            return self.__cplx_roots
        elif self.__cplx_roots_prec < prec:
            old_rts = self.__cplx_roots
            new_crv = self.complex_genus_2_curve(prec) #increases complex precision
            CC = self.complex_field()
            new_rts = Sequence(new_crv.roots(multiplicities = False), universe = CC)
            self.__cplx_roots = [0 for i in range(7)]
            for i in range(7):
                nb = 0
                for j in range(7):
                    if new_rts[i].overlaps(old_rts[j]):
                        nb++
                        self.__cplx_roots[j] = new_rts[i]
                if nb != 1:
                    raise ValueError("Zero or several overlaps between Weierstrass points when increasing precision from {} to {}".format(self.__cplx_roots_prec, prec))
            self.__cplx_roots_prec = prec
            self.__cplx_roots = Sequence(self.__cplx_roots, universe = CC)
            return self.__cplx_roots
        else:
            CC = ComplexBallField(prec)
            return Sequence(self.__cplx_roots, universe = CC)

    @cached
    def thomae_signs(self):
        #call C library. Return a permutation of the roots and a set of sign choices

    def rosenhain_invariants(self, prec):
        rts = self.complex_weierstrass_points(prec)
        #call C library with perm

    def theta4(self, prec):
        if self.is_geometrically_split():
            pass
        else:
            ros = self.rosenhain_invariants(prec);
            #call C library

    def theta2(self, prec):
        if self.is_geometrically_split():
            pass
        else:
            ros = self.rosenhain_invariants(prec)
            th4 = self.theta4(prec)
            #call C library with signs

    def complex_j_invariants(self, prec):
        CC = ComplexBallField(prec)
        pol = self.rational_invariants().j_invariants_quadratic_equation()
        delta = pol.discriminant()
        if delta == 0:
            j = pol.roots(multiplicities = False)[0]
            return Sequence([j, j], universe = CC)
        else:
            j = pol.roots(ComplexBallField(prec), multiplicities = False)
            return Sequence(j, universe = CC)

    def reduced_small_period_matrix(self, prec):
        if self.__cplx_periods is None or self.__cplx_periods_prec < prec:
            if self.is_geometrically_split():
                j1, j2 = self.complex_j_invariants(prec)
                tau1 = #call C library
                tau2 = #call C library
                CC = j1.parent()
                self.__cplx_periods_prec = prec
                self.__cplx_periods = Matrix(CC, [[tau1, 0], [0, tau2]])
            else:
                th2 = self.theta2(prec)
                #call C library

    

