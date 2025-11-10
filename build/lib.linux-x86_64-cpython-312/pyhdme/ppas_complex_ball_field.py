r"""
Principally polarized abelian surfaces over QQ and complex period matrices

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

from sage.structure.sage_object import SageObject
from sage.rings.complex_arb import ComplexBallField
from pyhdme.ppas_invariants import PPASInvariants

# Copyright 2025 Jean Kieffer
# See LICENSE file for license details.

class PPASComplexBallField(SageObject):
    r"""b

    Create a data structure encoding a principally polarized abelian surface
    over a complex ball field.
    """

    def P1_simple_roots(pol, prec = None):
        if pol.degree() > 6:
            raise ValueError("Polynomial must have degree at most 6")
        if prec is None:
            try:
                prec = pol.parent().base_ring().precision()
            except AttributeError:
                raise ValueError("Precision must be specified if not supplied by the base ring")
        CC = ComplexBallField(prec)
        if pol.coefficient(6) == 0:
            assert pol.degree() == 5, "Polynomial must have degree 5 or 6"
            rts = pol.roots()
            return [(x, CC(1)) for x in rts] + [(CC(1), CC(0))]
        elif pol.coefficient(6).accuracy() > 8/10 * prec:
            rts = pol.roots()
            return [(x, CC(1)) for x in rts]
        else:
            # Find an integer 0 <= i < 10 with P(i) nonzero
            i0 = -1
            for i in range(10):
                if pol(i).accuracy() > 8/10 * prec:
                    i0 = i
                    break
            if i0 == -1:
                raise ValueError("Could not find a nonzero value: polynomial too imprecise?")
            # Make transformation to ensure pol has degree 6, ensuring pol'(oo) = pol(i0)
            x = pol.parent().gen()
            pol = pol(i0 * x - 1).reverse(6)
            rts = pol.roots()
            return [(i0 * x - 1, x) for x in rts]

    def P1_no_overlap(x, y):
        x1, x2 = x
        y1, y2 = y
        t = x1 * y2 - y1 * x2
        if t.parent().is_exact():
            return t != 0
        else:
            return t.is_nonzero()

    def __init__(self, data, prec = None, low_precision = None):

        self.__prec = None
        self.__invs = None
        self.__crv = None
        self.__weierstrass = None
        self.__weierstrass_low = None
        self.__weierstrass_order = None
        self.__theta2_signs = None
        self.__theta2 = None
        self.__periods = None
        self.__gl2 = None
        self.__gl2_det = None

        # Initialize from period matrix
        if isinstance(data, Matrix):
            if prec is None:
                prec = data.base_ring().precision()
            CC = ComplexBallField(prec)
            data = Matrix(CC, 2, 2, data)
            self.__prec = prec
            self.__periods = data

        # Initialize from PPASInvariants structure
        elif isinstance(data, PPASInvariants):
            #todo: check automorphism group has been computed
            if prec is None:
                prec = data.base_ring().precision()
            CC = ComplexBallField(prec)
            data = data.change_ring(CC)
            self.__prec = prec
            self.__invs = data
            if isinstance(low_precision, PPASComplexBallField):
                #todo: compatibility checks
                self.__weierstrass_ordering = low_precision.__weierstrass_ordering
                self.__theta2_signs = low_precision.__theta2_signs
                self.__weierstrass_low = low_precision.__weierstrass

    def precision(self):
        return self.__prec

    def base_ring(self):
        return ComplexBallField(self.prec)

    def invariants(self):
        if self.__invs is None:
            self.__theta2 = 0 #call C library
            modular = 0 #call C library
            self.__invs = PPASInvariants(modular)
            self.__gl2_det = self.base_ring()(1)
        return self.__invs

    def has_generic_automorphisms(self):
        try:
            return self.invariants().geometric_automorphism_group_order() == 2
        except ValueError:
            #call C library
            pass

    def genus_2_curve_equation(self):
        if self.__crv is None:
            if self.__periods is None:
                self.__crv = self.invariants().genus_2_curve_equation() #this will fail if automorphisms were not computed
            else:
                self.__crv = 0 #call C library
                self.__gl2 = Matrix(self.base_ring(), 2, 2, 1)
                self.__gl2_det = self.base_ring()(1)
        return self.__crv

    def weierstrass_points(self):
        if self.__weierstrass is None:
            crv = self.genus_2_curve_equation()
            rts = crv.roots(multiplicities = False)
            if len(rts) != 6:
                raise ValueError("Unable to isolate six roots") #what do we do?
            if self.__weierstrass_low is None:
                # Check the roots do not overlap
                for i in range(7):
                    for j in range(i + 1, 7):
                        if rts[i].overlaps(rts[j]):
                            raise ValueError("Overlapping roots")
                self.__weierstrass = rts
            else:
                # Check the roots overlap exactly one low-precision root, and
                # reorder them
                self.__weierstrass = [0 for i in range(7)]
                for i in range(7):
                    nb = 0
                    for j in range(7):
                        if rts[i].overlaps(self.__weierstrass_low[j]):
                            nb += 1
                            self.__weierstrass[j] = rst[i]
                    if nb != 1:
                        raise ValueError("Zero or several overlaps between Weierstrass points")
        return self.__weierstrass

    def thomae_signs(self):
        if self.__weierstrass_order is None:
            # call C library
            pass
        return self.__weierstrass_order, self.__theta2_signs

    def rosenhain_invariants(self):
        rts = self.weierstrass_points()
        perm = self.thomae_signs()[0]
        # call C library

    def theta4(self):
        if self.is_geometrically_split():
            pass
        else:
            ros = self.rosenhain_invariants(prec);
            #call C library

    def theta2(self):
        if self.is_geometrically_split():
            pass
        else:
            ros = self.rosenhain_invariants(prec)
            th4 = self.theta4(prec)
            #call C library with signs

    def reduced_small_period_matrix(self):
        if self.__cplx_periods is None:
            if self.is_geometrically_split():
                j1, j2 = self.complex_j_invariants(prec)
                tau1 = 0 #call C library
                tau2 = 0 #call C library
                CC = j1.parent()
                self.__cplx_periods_prec = prec
                self.__cplx_periods = Matrix(CC, [[tau1, 0], [0, tau2]])
            else:
                th2 = self.theta2(prec)
                #call C library
        pass

    def GL2_on_differentials(self):
        if self.__gl2 is None:
            crv1 = 0 #call C library
            crv2 = self.genus_2_curve_equation()
            self.__gl2 = 0 #find isomorphism
            self.__gl2_det = self.__gl2.determinant()
        return self.__gl2

    def det_on_differentials(self):
        if self.__gl2_det is None:
            self.__gl2_det = 0 #find rescaling
        return self.__gl2_det
