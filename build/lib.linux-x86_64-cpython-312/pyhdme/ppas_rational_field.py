r"""
Principally polarized abelian surfaces over QQ and complex period matrices

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

from sage.rings.complex_arb import ComplexBallField
from pyhdme.ppas_invariants import PPASInvariants
from pyhdme.ppax_complex_ball_field import PPASComplexBallField

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
        self.__cplx_ppas = PPASComplexBallField(self.__cplx_invariants, prec = self.__cplx_prec)

    def rational_invariants(self):
        return self.__qq_invariants

    def increase_prec(self, prec):
        if self.__cplx_prec < prec:
            self.__cplx_prec = prec
            # this includes base-changing the curve equation if it was already computed.
            self.__cplx_invariants = self.rational_invariants().change_ring(CC)
            # this keeps thomae signs, etc. if they were already computed.
            self.__cplx_ppas = PPASComplexBallField(self.__cplx_invariants, prec = self.__cplx_prec,
                                                    low_precision = self.__cplx_ppas)
    def complex_invariants(self, prec = self.__cplx_prec):
        if self.__cplx_prec < prec:
            self.increase_prec(prec)
        elif self.__cplx_prec > prec:
            CC = ComplexBallField(prec)
            return self.__cplx_invariants.change_ring(CC)
        return self.__cplx_invariant

    def complex_ppas(self, prec = self.__cplx_prec):
        if self.__cplx_prec < prec:
            self.increase_prec(prec)
        elif self.__cplx_prec > prec:
            CC = ComplexBallField(prec)
            return self.__cplx_ppas.change_ring(CC)
        return self.__cplx_ppas

    def is_geometrically_split(self):
        return self.rational_invariants().is_geometrically_split()

    def complex_field(self):
        return self.complex_invariants().base_ring()

    def complex_genus_2_curve(self, prec = self.__cplx_prec):
        if self.is_geometrically_split():
            raise ValueError("This PPAS is geometrically split")
        return self.complex_invariants(prec).genus_2_curve_equation()

    def complex_weierstrass_points(self, prec = self.__cplx_prec):
        p = max(prec, 32)
        while True:
            try:
                return self.complex_ppas(p).weierstrass_points()
            except ValueError:
                p *= 2

    def reduced_small_period_matrix(self, prec = self.__cplx_prec):
        p = max(prec, 32)
        while True:
            try:
                return self.complex_ppas(p).reduced_small_period_matrix()
            except ValueError:
                p *= 2

    def rosenhain_invariants(self, prec = self.__cplx_prec):
        self.reduced_small_period_matrix()
        return self.complex_ppas().rosenhain_invariants()

    def theta2(self, prec = self.__cplx_prec):
        self.reduced_small_period_matrix()
        return self.complex_ppas().theta2()

    def complex_j_invariants(self, prec = self.__cplx_prec):
        CC = ComplexBallField(prec)
        pol = self.rational_invariants().j_invariants_quadratic_equation()
        delta = pol.discriminant()
        if delta == 0:
            j = pol.roots(multiplicities = False)[0]
            return Sequence([j, j], universe = CC)
        else:
            j = pol.roots(ComplexBallField(prec), multiplicities = False)
            return Sequence(j, universe = CC)

