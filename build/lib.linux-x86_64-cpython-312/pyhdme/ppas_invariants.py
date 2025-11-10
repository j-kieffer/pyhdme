r"""
Invariants and reconstruction of principally polarized abelian surfaces

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

# Copyright 2025 Jean Kieffer
# See LICENSE file for license details.

from sage.misc.cachefunc import cached_method
from sage.misc.functional import sqrt
from sage.arith.misc import gcd, xgcd
from sage.arith.functions import lcm
from sage.structure.sage_object import SageObject
from sage.structure.sequence import Sequence
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.rings.polynomial.polynomial_element import Polynomial
from sage.schemes.plane_conics.constructor import Conic
from sage.schemes.elliptic_curves.constructor import EllipticCurve, EllipticCurve_from_j
from sage.schemes.elliptic_curves.ell_generic import EllipticCurve_generic

class PPASInvariants(SageObject):
    r"""
    Create a data structure allowing conversions between different types of
    invariants of principally polarized abelian surfaces over a field.

    We require that 2, 3 and 5 are invertible in the base field.

    The different sets of invariants are as follows:

    - The "modular invariants", which correspond to the following Siegel
      modular forms with integral Fourier expansions:

      \psi_4 = 1 + 240 * (q_1 + q_2) + ...

      \psi_6 = 1 - 504 * (q_1 + q_2) + ...

      \chi_{10} = (q_3 - 2 + q_3^{-1}) + ...

      \chi_{12} = (q_3 + 10 + q_3^{-1}) + ...

      These invariants make sense both for Jacobians of genus 2 curves and for
      products of elliptic curves.

    - The classical Igusa--Clebsch invariants I_2, I_4, I_6, I_{10}.

    - The classical Clebsch invariants A, B, C, D.

    - The modified Igusa--Clebsch invariants I_4, I_6', I_{10}, I_{12}.

    - The absolute Igusa invariants, defined as follows:

      j_1 = I_4 * I_6' / I_{10},

      j_2 = I_4^2 * I_{12} / I_{10}^2,

      j_3 = I_4^5 / I_{10}^2.

    Some properties of principally polarized abelian surfaces that are easily
    obtained from invariants (such as automorphism groups) are also available.
    A :class:`ValueError` is raised when the conversion doesn't make sense
    (e.g. when asking for the Igusa invariants of a product of elliptic
    curves).

    """

    def ambient_field(F):
        r"""
        Return a minimal base field F' that contains the ring F, and raises
        a :class:`ValueError` if F' has characteristic 2, 3, or 5.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants.ambient_field(ZZ)
            Rational Field
            sage: PPASInvariants.ambient_field(FiniteField(5))
            Traceback (most recent call last):
            ...
            ValueError: Base field cannot have characteristic 2, 3 or 5

        """
        if not F.is_field():
            F = F.fraction_field()
        if F(2) == 0 or F(3) == 0 or F(5) == 0:
            raise ValueError("Base field cannot have characteristic 2, 3 or 5")
        return F

    def modular_from_modified_igusa(vec):
        r"""
        Return the modular invariants associated to the given modified Igusa--Clebsch invariants.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "ModifiedIgusaClebsch").modular_invariants()
            sage: vec == PPASInvariants.modular_from_modified_igusa([1, 2, 3, 4])
            True

        """
        I4, I6p, I10, I12 = vec
        return [I4 / 4, I6p / 4, -I10 / 2**12, I12 / 2**15]

    def modified_igusa_clebsch_from_igusa_clebsch(vec):
        r"""
        Return the modified Igusa--Clebsch invariants associated to the given Igusa--Clebsch invariants.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "IgusaClebsch").modified_igusa_clebsch_invariants()
            sage: vec == PPASInvariants.modified_igusa_clebsch_from_igusa_clebsch([1, 2, 3, 4])
            True

        """
        I2, I4, I6, I10 = vec
        return [I4, (I2 * I4 - 3 * I6) / 2, I10, I2 * I10]

    def igusa_clebsch_from_clebsch(vec):
        r"""
        Return the Igusa--Clebsch invariants associated to the given Clebsch invariants.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "Clebsch").igusa_clebsch_invariants()
            sage: vec == PPASInvariants.igusa_clebsch_from_clebsch([1, 2, 3, 4])
            True

        """
        A, B, C, D = vec
        I2 = - 120 * A
        I4 = - 720 * A**2 + 6750 * B
        I6 = 8640 * A**3 - 108000 * A * B + 202500 * C
        I10 = - 62208 * A**5 + 972000 * A**3 * B + 1620000 * A**2 * C - 3037500 * A * B**2 - 6075000 * B * C - 4556250 * D
        return [I2, I4, I6, I10]

    def modified_igusa_clebsch_from_curve(curve):
        r"""
        Return the modified Igusa--Clebsch invariants from the given vector
        of curve coefficients. This should be rewritten in terms of
        transvectants of binary forms for efficiency.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: R.<t> = PolynomialRing(QQ)
            sage: vec = PPASInvariants(t^6 - 4*t^4 + 7).modified_igusa_clebsch_invariants()
            sage: vec == PPASInvariants.modified_igusa_clebsch_from_curve(t^6 - 4*t^4 + 7)
            True

        """
        a, b, c, d, e, f, g = [curve.coefficient(6 - i) for i in range(7)]
        I2 = -240 * g * a + (40 * f * b + (-16 * e * c + 6 * d**2))
        I4 = 1620 * g**2 * a**2 + (-540 * g * f * b + ((-504 * g * e + 300 * f**2) * c + (324 * g * d**2 - 180 * f * e * d + 48 * e**3))) * a + ((300 * g * e - 80 * f**2) * b**2 + ((-180 * g * d + 4 * f * e) * c + (36 * f * d**2 - 12 * e**2 * d)) * b + (48 * g * c**3 + (-12 * f * d + 4 * e**2) * c**2))
        I6p = -14580 * g**3 * a**3 + (7290 * g**2 * f * b + ((16524 * g**2 * e -8100 * g * f**2) * c + (-18954 * g**2 * d**2 + (17010 * g * f * e - 3375 * f**3) * d+ (-5616 * g * e**3 + 1350 * f**2 * e**2)))) * a**2 + ((-8100 * g**2 * e +2160 * g * f**2) * b**2 + ((17010 * g**2 * d + (-11448 * g * f * e +3600 * f**3)) * c + (-2187 * g * f * d**2 + (2754 * g * e**2 - 810 * f**2 * e) * d +36 * f * e**3)) * b + (-5616 * g**2 * c**3 + (2754 * g * f * d + (2916 * g * e**2 -1440 * f**2 * e)) * c**2 + ((-3402 * g * e + 405 * f**2) * d**2 + 702 * f * e**2 * d- 144 * e**4) * c + (729 * g * d**4 - 243 * f * e * d**3 + 54 * e**3 * d**2))) * a +((-3375 * g**2 * d + (3600 * g * f * e - 1120 * f**3)) * b**3 + (1350 * g**2 * c**2+ (-810 * g * f * d + (-1440 * g * e**2 + 624 * f**2 * e)) * c + ((405 * g * e +216 * f**2) * d**2 - 279 * f * e**2 * d + 54 * e**4)) * b**2 + (36 * g * f * c**3 +((702 * g * e - 279 * f**2) * d + 6 * f * e**2) * c**2 + (-243 * g * d**3 +81 * f * e * d**2 - 18 * e**3 * d) * c) * b + ((-144 * g * e + 54 * f**2) * c**4 +(54 * g * d**2 - 18 * f * e * d + 4 * e**3) * c**3))
        I10 = -46656 * g**5 * a**5 + (38880 * g**4 * f * b + ((62208 * g**4 * e -32400 * g**3 * f**2) * c + (34992 * g**4 * d**2 + (-77760 * g**3 * f * e +27000 * g**2 * f**3) * d + (-13824 * g**3 * e**3 + 43200 * g**2 * f**2 * e**2 -22500 * g * f**4 * e + 3125 * f**6)))) * a**4 + ((-32400 * g**4 * e +540 * g**3 * f**2) * b**2 + ((-77760 * g**4 * d + (31968 * g**3 * f * e -1800 * g**2 * f**3)) * c + (15552 * g**3 * f * d**2 + (46656 * g**3 * e**2 -31320 * g**2 * f**2 * e + 2250 * g * f**4) * d + (-21888 * g**2 * f * e**3 +15600 * g * f**3 * e**2 - 2500 * f**5 * e))) * b + (-13824 * g**4 * c**3 +(46656 * g**3 * f * d + (-17280 * g**3 * e**2 - 6480 * g**2 * f**2 * e +1500 * g * f**4)) * c**2 + ((3888 * g**3 * e - 27540 * g**2 * f**2) * d**2 +(-3456 * g**2 * f * e**2 + 19800 * g * f**3 * e - 3750 * f**5) * d +(9216 * g**2 * e**4 - 10560 * g * f**2 * e**3 + 2000 * f**4 * e**2)) * c +(-8748 * g**3 * d**4 + (21384 * g**2 * f * e - 1350 * g * f**3) * d**3 +(-8640 * g**2 * e**3 - 9720 * g * f**2 * e**2 + 2250 * f**4 * e) * d**2 +(6912 * g * f * e**4 - 1600 * f**3 * e**3) * d + (-1024 * g * e**6 +256 * f**2 * e**5)))) * a**3 + ((27000 * g**4 * d + (-1800 * g**3 * f * e +410 * g**2 * f**3)) * b**3 + (43200 * g**4 * c**2 + (-31320 * g**3 * f * d +(-6480 * g**3 * e**2 + 8748 * g**2 * f**2 * e - 1700 * g * f**4)) * c +((-27540 * g**3 * e + 15417 * g**2 * f**2) * d**2 + (16632 * g**2 * f * e**2 -12330 * g * f**3 * e + 2000 * f**5) * d + (-192 * g**2 * e**4 + 248 * g * f**2 * e**3- 50 * f**4 * e**2))) * b**2 + (-21888 * g**3 * f * c**3 + ((-3456 * g**3 * e +16632 * g**2 * f**2) * d + (15264 * g**2 * f * e**2 - 13040 * g * f**3 * e +2250 * f**5)) * c**2 + (21384 * g**3 * d**3 + (-22896 * g**2 * f * e +1980 * g * f**3) * d**2 + (-5760 * g**2 * e**3 + 10152 * g * f**2 * e**2 -2050 * f**4 * e) * d + (-640 * g * f * e**4 + 160 * f**3 * e**3)) * c +(-6318 * g**2 * f * d**4 + (5832 * g**2 * e**2 + 3942 * g * f**2 * e -900 * f**4) * d**3 + (-4464 * g * f * e**3 + 1020 * f**3 * e**2) * d**2 +(768 * g * e**5 - 192 * f**2 * e**4) * d)) * b + ((9216 * g**3 * e -192 * g**2 * f**2) * c**4 + (-8640 * g**3 * d**2 + (-5760 * g**2 * f * e -120 * g * f**3) * d + (-4352 * g**2 * e**3 + 4816 * g * f**2 * e**2 -900 * f**4 * e)) * c**3 + (5832 * g**2 * f * d**3 + (8208 * g**2 * e**2 -4536 * g * f**2 * e + 825 * f**4) * d**2 + (-2496 * g * f * e**3 +560 * f**3 * e**2) * d + (512 * g * e**5 - 128 * f**2 * e**4)) * c**2 +((-4860 * g**2 * e + 162 * g * f**2) * d**4 + (2808 * g * f * e**2 -630 * f**3 * e) * d**3 + (-576 * g * e**4 + 144 * f**2 * e**3) * d**2) * c +(729 * g**2 * d**6 + (-486 * g * f * e + 108 * f**3) * d**5 + (108 * g * e**3 -27 * f**2 * e**2) * d**4))) * a**2 + ((-22500 * g**4 * c + (2250 * g**3 * f * d +(1500 * g**3 * e**2 - 1700 * g**2 * f**2 * e + 320 * g * f**4))) * b**4 +(15600 * g**3 * f * c**2 + ((19800 * g**3 * e - 12330 * g**2 * f**2) * d +(-13040 * g**2 * f * e**2 + 9768 * g * f**3 * e - 1600 * f**5)) * c +(-1350 * g**3 * d**3 + (1980 * g**2 * f * e - 208 * g * f**3) * d**2 +(-120 * g**2 * e**3 - 682 * g * f**2 * e**2 + 160 * f**4 * e) * d + (144 * g * f * e**4- 36 * f**3 * e**3))) * b**3 + ((-10560 * g**3 * e + 248 * g**2 * f**2) * c**3 +(-9720 * g**3 * d**2 + (10152 * g**2 * f * e - 682 * g * f**3) * d +(4816 * g**2 * e**3 - 5428 * g * f**2 * e**2 + 1020 * f**4 * e)) * c**2 +(3942 * g**2 * f * d**3 + (-4536 * g**2 * e**2 - 2412 * g * f**2 * e +560 * f**4) * d**2 + (3272 * g * f * e**3 - 746 * f**3 * e**2) * d + (-576 * g * e**5+ 144 * f**2 * e**4)) * c + (162 * g**2 * e * d**4 + (-108 * g * f * e**2 +24 * f**3 * e) * d**3 + (24 * g * e**4 - 6 * f**2 * e**3) * d**2)) * b**2 +((6912 * g**3 * d + (-640 * g**2 * f * e + 144 * g * f**3)) * c**4 +(-4464 * g**2 * f * d**2 + (-2496 * g**2 * e**2 + 3272 * g * f**2 * e -630 * f**4) * d + (-96 * g * f * e**3 + 24 * f**3 * e**2)) * c**3 + ((2808 * g**2 * e- 108 * g * f**2) * d**3 + (-1584 * g * f * e**2 + 356 * f**3 * e) * d**2 +(320 * g * e**4 - 80 * f**2 * e**3) * d) * c**2 + (-486 * g**2 * d**5 +(324 * g * f * e - 72 * f**3) * d**4 + (-72 * g * e**3 +18 * f**2 * e**2) * d**3) * c) * b + (-1024 * g**3 * c**6 + (768 * g**2 * f * d +(512 * g**2 * e**2 - 576 * g * f**2 * e + 108 * f**4)) * c**5 + ((-576 * g**2 * e +24 * g * f**2) * d**2 + (320 * g * f * e**2 - 72 * f**3 * e) * d + (-64 * g * e**4 +16 * f**2 * e**3)) * c**4 + (108 * g**2 * d**4 + (-72 * g * f * e + 16 * f**3) * d**3 +(16 * g * e**3 - 4 * f**2 * e**2) * d**2) * c**3)) * a + (3125 * g**4 * b**6 +(-2500 * g**3 * f * c + ((-3750 * g**3 * e + 2000 * g**2 * f**2) * d +(2250 * g**2 * f * e**2 - 1600 * g * f**3 * e + 256 * f**5))) * b**5 +((2000 * g**3 * e - 50 * g**2 * f**2) * c**2 + (2250 * g**3 * d**2 +(-2050 * g**2 * f * e + 160 * g * f**3) * d + (-900 * g**2 * e**3 +1020 * g * f**2 * e**2 - 192 * f**4 * e)) * c + (-900 * g**2 * f * d**3 +(825 * g**2 * e**2 + 560 * g * f**2 * e - 128 * f**4) * d**2 + (-630 * g * f * e**3 +144 * f**3 * e**2) * d + (108 * g * e**5 - 27 * f**2 * e**4))) * b**4 +((-1600 * g**3 * d + (160 * g**2 * f * e - 36 * g * f**3)) * c**3 +(1020 * g**2 * f * d**2 + (560 * g**2 * e**2 - 746 * g * f**2 * e + 144 * f**4) * d +(24 * g * f * e**3 - 6 * f**3 * e**2)) * c**2 + ((-630 * g**2 * e +24 * g * f**2) * d**3 + (356 * g * f * e**2 - 80 * f**3 * e) * d**2 + (-72 * g * e**4 +18 * f**2 * e**3) * d) * c + (108 * g**2 * d**5 + (-72 * g * f * e + 16 * f**3) * d**4+ (16 * g * e**3 - 4 * f**2 * e**2) * d**3)) * b**3 + (256 * g**3 * c**5 +(-192 * g**2 * f * d + (-128 * g**2 * e**2 + 144 * g * f**2 * e - 27 * f**4)) * c**4 +((144 * g**2 * e - 6 * g * f**2) * d**2 + (-80 * g * f * e**2 + 18 * f**3 * e) * d +(16 * g * e**4 - 4 * f**2 * e**3)) * c**3 + (-27 * g**2 * d**4 + (18 * g * f * e -4 * f**3) * d**3 + (-4 * g * e**3 + f**2 * e**2) * d**2) * c**2) * b**2)
        return [I4, I6p, I10, I2 * I10]

    def parametrize_conic(pt, conic, t):
        r"""
        Return a parametrization of the given conic from the given base point,
        using t as variable.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: conic = [1, 1, -1, 0, 0, 0]                # x^2 + y^2 - z^2
            sage: pt = [0, -1, 1]
            sage: R.<t> = PolynomialRing(QQ)
            sage: x, y, z = PPASInvariants.parametrize_conic(pt, conic, t)
            sage: x^2 + y^2 - z^2 == 0
            True
            sage: x.degree() > 0
            True

        """
        x0, y0, z0 = pt
        c11, c22, c33, c23, c31, c12 = conic

        #Enforce x0 != 0
        if x0 == 0:
            if y0 != 0:
                y, z, x = PPASInvariants.parametrize_conic([y0, z0, x0],
                                                           [c22, c33, c11, c31, c12, c23], t)
                return [x, y, z]
            elif z0 != 0:
                z, x, y = PPASInvariants.parametrize_conic([z0, x0, y0],
                                                           [c33, c11, c22, c12, c23, c13], t)
                return [x, y, z]
            else:
                raise ValueError("Conic point does not have any nonzero coordinates")

        R = PolynomialRing(t.parent(), "u")
        u = R.gen()
        x = x0
        y = y0 + u * t
        z = z0 + u
        substitution = c11 * x**2 + c22 * y**2 + c33 * z**2 + 2 * c23 * y * z + 2 * c31 * x * z + 2 * c12 * x * y
        a = substitution.coefficient(2) # in u
        b = substitution.coefficient(1) # in u
        return [x0 * a, y0 * a - t * b, z0 * a - b]

    def find_rescaling(R, v1, v2, e, w):
        r"""
        Given two vectors v1 and v2, find an element lambda in the ring R such
        that \prod_i v1_i^{e_i} = lambda^k * \prod_i v2_i^{e_i}, where we set
        k = \sum_i e_i w_i. We require k > 0.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: v1 = [0, 1, 1]
            sage: v2 = [0, -1, 1]
            sage: e = [0, -1, 1]
            sage: w = [1, 3, 6]
            sage: PPASInvariants.find_rescaling(Integers(), v1, v2, e, w) == -1
            True

        """
        a1 = 1
        a2 = 1
        k = 0
        n = len(v1)
        if n != len(v2) or n != len(e) or n != len(w):
            raise ValueError("Input vectors must be of the same length")
        for i in range(n):
            a1 *= v1[i] ** e[i]
            a2 *= v2[i] ** e[i]
            k += e[i] * w[i]
        if k <= 0:
            raise ValueError("Minimal weight must be positive")
        x = a1 / a2
        try:
            x = x.nth_root(k)
        except (TypeError, ValueError):
            raise ValueError("Could not extract a {}th root of {} in {}".format(k, x, R))
        return x

    def __init__(self, data, inv_type = "Modular"):
        r"""
        Initialize a PPASInvariants data structure. The input can be one
        of the following:

        - A polynomial f of degree 5 or 6 encoding the genus 2 curve y^2 = f(x),

        - A tuple of 7 coefficients a_6, ..., a_0, encoding the polynomial
          f = a_6 x^6 + ... + a_0 as in the first item,

        - A pair of elliptic curves,

        - A triple of absolute Igusa invariants, or

        - A tuple of 4 invariants, whose type is specified by inv_type. The
          possible values are: "Modular" (default), "IgusaClebsch",
          "ModifiedIgusaClebsch", and "Clebsch".

        The corresponding modular invariants are then computed.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: R.<x> = PolynomialRing(QQ)
            sage: PPASInvariants(x^6 + 2*x^4 - x^2 + 7).modular_invariants()
            [22273, -2171249, 22453767/64, 2312738001/32]
            sage: PPASInvariants([1, 0, 2, 0, -1, 0, 7]).modular_invariants()
            [22273, -2171249, 22453767/64, 2312738001/32]
            sage: PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).modular_invariants()
            [6912, 5971968, 0, 185794560]
            sage: PPASInvariants([1, 2, 3]).modular_invariants()
            [3/4, 3/4, -9/4096, 9/16384]
            sage: PPASInvariants([1, 2, 3, 4]).modular_invariants()
            [1, 2, 3, 4]
            sage: PPASInvariants([1, 2, 3, 4], "IgusaClebsch").modular_invariants()
            [1/2, -7/8, -1/1024, 1/8192]
            sage: PPASInvariants([1, 2, 3, 4], "ModifiedIgusaClebsch").modular_invariants()
            [1/4, 1/2, -3/4096, 1/8192]
            sage: PPASInvariants([1, 2, 3, 4], "Clebsch").modular_invariants()
            [3195, -683505/2, 7510401/512, 112656015/512]

        """

        self.__g2_curve = None
        self.__ell_curves = None

        self.__modular = None
        self.__ic = None
        self.__clebsch = None
        self.__ic_mod = None
        self.__abs_igusa = None

        self.__mestre_U = None
        self.__mestre_line = None
        self.__aut_gp_order = None
        self.__bolza_a2 = None
        self.__min_wt = None

        if isinstance(data, Polynomial):
            F = data.parent().base_ring()
            F = PPASInvariants.ambient_field(F)
            self.__g2_curve = data.base_extend(F)
            self.__ic_mod = Sequence(PPASInvariants.modified_igusa_clebsch_from_curve(self.__g2_curve),
                                     universe = F)
            self.__modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.__ic_mod),
                                      universe = F)

        elif isinstance(data, list):

            if len(data) == 2:
                E1 = data[0]
                E2 = data[1]
                if not (isinstance(E1, EllipticCurve_generic) and isinstance(E2, EllipticCurve_generic)):
                    raise TypeError("Input must be a pair of elliptic curves")
                F = Sequence(E1.a_invariants() + E2.a_invariants()).universe()
                F = PPASInvariants.ambient_field(F)
                self.__ell_curves = [E1.change_ring(F), E2.change_ring(F)]
                c41, c61 = E1.c_invariants()
                c42, c62 = E2.c_invariants()
                delta1 = E1.discriminant()
                delta2 = E2.discriminant()
                self.__modular = Sequence([c41 * c42, c61 * c62, 0, 12 * delta1 * delta2], universe = F)

            elif len(data) == 3:
                F = Sequence(data).universe()
                F = PPASInvariants.ambient_field(F)
                self.__abs_igusa = Sequence(data, universe = F)
                self.__ic_mod = Sequence([data[2], data[0] * data[2], data[2]**2,
                                          data[1] * data[2]**2], universe = F)
                self.__modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.__ic_mod),
                                          universe = F)

            elif len(data) == 4:
                F = Sequence(data).universe()
                F = PPASInvariants.ambient_field(F)
                data = Sequence(data, universe = F)

                if inv_type == "Clebsch":
                    self.__clebsch = data
                    data = Sequence(PPASInvariants.igusa_clebsch_from_clebsch(self.__clebsch),
                                    universe = F)
                if inv_type in ["Clebsch", "IgusaClebsch"]:
                    self.__ic = data
                    data = Sequence(PPASInvariants.modified_igusa_clebsch_from_igusa_clebsch(self.__ic),
                                    universe = F)
                if inv_type in ["Clebsch", "IgusaClebsch", "ModifiedIgusaClebsch"]:
                    self.__ic_mod = data
                    self.__modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.__ic_mod),
                                            universe = F)
                elif inv_type == "Modular":
                    self.__modular = data
                else:
                    raise TypeError("Unknown invariant type: {}".format(inv_type))

            elif len(data) == 7:
                F = Sequence(data).universe()
                F = PPASInvariants.ambient_field(F)
                R = PolynomialRing(F, "x")
                self.__g2_curve = R(data).reverse(6)
                self.__ic_mod = Sequence(PPASInvariants.modified_igusa_clebsch_from_curve(self.__g2_curve),
                                         universe = F)
                self.__modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.__ic_mod),
                                          universe = F)

            else:
                raise TypeError("Invalid input length {}".format(len(data)))


        else:
            raise TypeError("Input must be either a polynomial, a pair of elliptic curves, or a list of coefficients or invariants")

    def __repr__(self):
        r"""
        Return a string representation of self.

        EXAMPLES:

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 3, 4])
            Modular invariants of a principally polarized abelian surface with values [1, 2, 3, 4] over Rational Field

        """
        return "Modular invariants of a principally polarized abelian surface with values {} over {}".format(self.modular_invariants(), self.base_ring())

    def __eq__(self, other):
        return self.modular_invariants() == other.modular_invariants()

    def base_ring(self):
        r"""
        Return the base ring of the principally polarized abelian surface. This
        is always a field not of characteristic 2, 3, or 5.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 3, 4]).base_ring()
            Rational Field

        """
        return self.modular_invariants().universe()

    def modular_invariants(self):
        r"""
        Return the modular invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 3, 4]).modular_invariants()
            [1, 2, 3, 4]

        """
        return self.__modular

    def modified_igusa_clebsch_invariants(self):
        r"""
        Return the modified Igusa--Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "ModifiedIgusaClebsch").modular_invariants()
            sage: PPASInvariants(vec).modified_igusa_clebsch_invariants()
            [1, 2, 3, 4]

        """
        if self.__ic_mod is None:
            m4, m6, m10, m12 = self.modular_invariants()
            self.__ic_mod = Sequence([4 * m4, 4 * m6, - 2**12 * m10, 2**15 * m12], universe = self.base_ring())
        return self.__ic_mod

    def is_geometrically_split(self):
        r"""
        Return True iff the provided invariants correspond to a geometrically split PPAS.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).is_geometrically_split()
            True
            sage: PPASInvariants([1, 2, 3, 4, 5, 6, 7]).is_geometrically_split()
            False

        """
        return self.modular_invariants()[2] == 0

    def igusa_clebsch_invariants(self):
        r"""
        Return the Igusa--Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "IgusaClebsch").modular_invariants()
            sage: PPASInvariants(vec).igusa_clebsch_invariants()
            [1, 2, 3, 4]
            sage: PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).igusa_clebsch_invariants()
            Traceback (most recent call last):
            ...
            ValueError: Igusa--Clebsch, Clebsch, or R2 invariants are not defined for geometrically split surfaces

        """
        if self.__ic is None:
            if self.is_geometrically_split():
                raise ValueError("Igusa--Clebsch, Clebsch, or R2 invariants are not defined for geometrically split surfaces")
            I4, I6p, I10, I12 = self.modified_igusa_clebsch_invariants()
            I2 = I12 / I10
            I6 = (2 * I6p - I2 * I4) / (-3)
            self.__ic = Sequence([I2, I4, I6, I10], universe = self.base_ring())
        return self.__ic

    def clebsch_invariants(self):
        r"""
        Return the Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4], "Clebsch").modular_invariants()
            sage: PPASInvariants(vec).clebsch_invariants()
            [1, 2, 3, 4]
            sage: PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).clebsch_invariants()
            Traceback (most recent call last):
            ...
            ValueError: Igusa--Clebsch, Clebsch, or R2 invariants are not defined for geometrically split surfaces

        """
        if self.__clebsch is None:
            I2, I4, I6, I10 = self.igusa_clebsch_invariants()
            A = -I2 / 120
            B = (I4 + 720 * A**2) / 6750
            C = (I6 - 8640 * A**3 + 108000 * A * B) / 202500
            D = (I10 + 62208 * A**5 - 972000 * A**3 * B - 1620000 * A**2 * C + 3037500 * A * B**2 + 6075000 * B * C) / (-4556250)
            self.__clebsch = Sequence([A, B, C, D], universe = self.base_ring())
        return self.__clebsch

    @cached_method
    def R2_invariant(self):
        r"""
        Return the invariant R^2, of weight 30, of the principally polarized abelian surface.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 3, 4]).R2_invariant()
            -701177130860721779769344/3
            sage: PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).R2_invariant()
            Traceback (most recent call last):
            ...
            ValueError: Igusa--Clebsch, Clebsch, or R2 invariants are not defined for geometrically split surfaces

        """
        a, b, c, d = self.igusa_clebsch_invariants()
        R2 = 125971200000 * d**3 + 236196 * d**2 * a**5 + 19245600 * d**2 * a**3 * b - 507384000 * d**2 * a * b**2 - 972 * d * a**6 * b**2 - 77436 * d * a**4 * b**3 + 592272 * d * a**2 * b**4 + a**7 * b**4 - 41472 * d * b**5 + 78 * a**5 * b**5 - 159 * a**3 * b**6 + 80 * a * b**7 - 104976000 * d**2 * a**2 * c + 2099520000 * d**2 * b * c + 5832 * d * a**5 * b * c + 870912 * d * a**3 * b**2 * c - 4743360 * d * a * b**3 * c - 12 * a**6 * b**3 * c - 1332 * a**4 * b**4 * c + 1728 * a**2 * b**5 * c - 384 * b**6 * c - 8748 * d * a**4 * c**2 - 3090960 * d * a**2 * b * c**2 + 9331200 * d * b**2 * c**2 + 54 * a**5 * b**2 * c**2 + 8910 * a**3 * b**3 * c**2 - 6048 * a * b**4 * c**2 + 3499200 * d * a * c**3 - 108 * a**4 * b * c**3 - 29376 * a**2 * b**2 * c**3 + 6912 * b**3 * c**3 + 81 * a**3 * c**4 + 47952 * a * b * c**4 - 31104 * c**5
        return R2

    @cached_method
    def igusa_integral_invariants(self):
        r"""
        Return the values of the 14 generators X_4, X_6, ..., X_{48} for the
        ring of even-weight Siegel modular forms over `\ZZ`, of weights 4, 6,
        10, 12, 12, 16, 18, 24, 28, 30, 36, 40, 42, 48 for the given PPAS.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1,2,3,4,5,6,7]).igusa_integral_invariants()
            sage: lcm([x.denominator() for x in vec]).divides(2**32)
            True

        """
        m4, m6, m10, m12 = self.modular_invariants()

        Y12 = (m4**3 - m6**2) / (2**6 * 3**3) + 2**4 * 3**2 * m12
        X16 = (m4 * m12 - m6 * m10) / 12
        X18 = (m6 * m12 - m4**2 * m10) / 12
        X24 = (m12**2 - m4 * m10**2) / 24
        X28 = (X24 * m4 - m10 * X18) / 6
        X30 = (m6 * X24 - X16 * m4 * m10) / 6
        X36 = (m12 * X24 - X16 * m10**2) / 18
        X40 = (m4 * X36 - m10 * X30) / 4
        X42 = (m12 * X30 - X28 * m4 * m10) / 12
        X48 = (m12 * X36 - X24**2) / 4
        return Sequence([m4, m6, m10, m12, Y12, X16, X18, X24, X28, X30, X36, X40, X42, X48],
                        universe = self.base_ring())

    def absolute_igusa_invariants(self):
        r"""
        Returns the absolute Igusa invariants j1, j2, j3 attached to the given
        PPASInvariants structure.

        EXAMPLES:

            sage: from pyhdme import PPASInvariants
            sage: v = PPASInvariants([1, 2, 3]).modular_invariants()
            sage: PPASInvariants(v).absolute_igusa_invariants()
            [1, 2, 3]

        """
        if self.__abs_igusa is None:
            if self.is_geometrically_split():
                raise ValueError("Absolute Igusa invariants are not defined for geometrically split surfaces")
            I4, I6p, I10, I12 = self.modified_igusa_clebsch_invariants()
            self.__abs_igusa = Sequence([I4 * I6p / I10, I4**2 * I12 / I10**2, I4**5 / I10**2],
                                        universe = self.base_ring())
        return self.__abs_igusa

    def minimal_weight_combination(self):
        r"""
        Returns a list of exponents [a, b, c, d] such that the product
        I_4^a I_6'^b I_{10}^c I_{12}^d is a nonvanishing monomial in the
        specified invariants of the smallest possible weight.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 3, 4]).minimal_weight_combination()
            [-1, 1, 0, 0]
            sage: PPASInvariants([0, 2, 3, 0]).minimal_weight_combination()
            [0, 2, -1, 0]
            sage: PPASInvariants([1, 0, 0, 4]).minimal_weight_combination()
            [1, 0, 0, 0]
            sage: PPASInvariants([0, 0, 3, 0]).minimal_weight_combination()
            [0, 0, 1, 0]

        """
        if not self.__min_wt is None:
            return self.__min_wt

        if not self.base_ring().is_exact():
            raise ValueError("Minimal weight combination is not implemented over inexact base fields")
        vec = self.modified_igusa_clebsch_invariants()
        weights = [4, 6, 10, 12]
        for i in range(4):
            if vec[i] == 0:
                weights[i] = 0
        res = xgcd(weights)
        self.__min_wt = list(res[1:len(res)])
        return self.__min_wt

    def _bolza_condition_19(self):
        r"""
        Returns True iff Bolza's condition 19 holds. If so, additionally return
        a^2 such that the given PPAS is isomorphic to the Jacobian of the genus
        2 curve y^2 = x^6 + a*x^3 + 1.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 5, 0, 0, 1])._bolza_condition_19()
            (True, 25)
            sage: PPASInvariants([1, 2, 3, 4, 5, 6, 7])._bolza_condition_19()
            (False, None)

        """
        A, B, C, D = self.clebsch_invariants()
        t1 = - B**3 + 6 * C**2
        t4 = A * B + 6 * C
        t2 = -2 * t4 * B + 9 * D
        t3 = -15 * C + 2 * A * B
        r = (t1 == 0) and (t2 == 0)
        if r:
            if (D == 0 or t3 == 0):
                raise ValueError("Unexpected vanishing")
            a2 = 100 * (-t4 / t3)
        else:
            a2 = None
        return r, a2

    def _bolza_condition_23(self):
        r"""
        Returns True iff Bolza's condition 23 holds. If so, additionally return
        a^2 such that the given PPAS is isomorphic to the Jacobian of the genus
        2 curve y^2 = x^5 + a*x^3 + x.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([0, 1, 0, 5, 0, 1, 0])._bolza_condition_23()
            (True, 25)
            sage: PPASInvariants([1, 2, 3, 4, 5, 6, 7])._bolza_condition_23()
            (False, None)

        """
        A, B, C, D = self.clebsch_invariants()
        t1 = 3 * B**2 * A - 6 * B * C + 4 * C * A**2 - 18 * D
        t2 = 4 * B**3 + 5 * C * B * A + 6 * C**2 - 3 * A * D
        t3 = 6 * C**2 - B**3
        r = (t1 == 0) and (t2 == 0)
        if r:
            if (D == 0 or t3 == 0):
                raise ValueError("Unexpected vanishing")
            a2 = 100 * (B**2 + A * C) / (2 * A**2 * B - 3 * A * C - 15 * B**2)
        else:
            a2 = None
        return r, a2

    def _bolza_a2(self):
        r"""
        Returns the square of the element a such that a genus 2 curve equation
        realizing the specified invariants is x^6 + a*x^3 + 1 or
        x^5 + a*x^3 + x. This requires the automorphism group to be of order 8
        or 12.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 7, 0, 0, 1])._bolza_a2()
            49
            sage: PPASInvariants([0, 1, 0, 7, 0, 1, 0])._bolza_a2()
            49

        """
        n = self.geometric_automorphism_group_order()
        if n != 8 and n != 12:
            raise ValueError("Automorphism group must have order 8 or 12")
        return self.__bolza_a2

    def geometric_automorphism_group_order(self):
        r"""
        Return the geometric automorphism group of the specified principally
        polarized abelian surface as an abstract group.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 2, 4, 4, 1]).geometric_automorphism_group_order()
            2
            sage: PPASInvariants([1, 0, 4, 2, 4, 0, 1]).geometric_automorphism_group_order()
            4
            sage: PPASInvariants([0, 8, -12, 4, 4, -4, 1]).geometric_automorphism_group_order()
            8
            sage: PPASInvariants([0, 4, 0, 0, 0, 0, 1]).geometric_automorphism_group_order()
            10
            sage: PPASInvariants([1, 4, 6, 2, 1, 2, 1]).geometric_automorphism_group_order()
            12
            sage: PPASInvariants([4, 0, 0, 0, 0, 0, 1]).geometric_automorphism_group_order()
            24
            sage: PPASInvariants([0, 1, 0, 0, 0, -1, 0]).geometric_automorphism_group_order()
            48

        """
        if not self.__aut_gp_order is None:
            return self.__aut_gp_order

        if not self.base_ring().is_exact():
            raise NotImplementedError("Automorphism group computation is not currently implemented over inexact fields")
        if self.is_geometrically_split():
            raise NotImplementedError("Automorphism groups are not currently implemented for products of elliptic curves")

        R2 = self.R2_invariant()
        A, B, C, D = self.clebsch_invariants()
        if R2 != 0 and (A != 0 or B != 0 or C != 0):
            self.__aut_gp_order = 2
        elif R2 != 0:
            self.__aut_gp_order = 10
        elif B == 0 and C == 0 and D == 0:
            self.__aut_gp_order = 48
        elif 6 * B - A**2 == 0 and 6 * C + A * B == 0 and D == 0:
            self.__aut_gp_order = 24
        else:
            r, a2 = self._bolza_condition_19()
            if r:
                self.__bolza_a2 = a2
                self.__aut_gp_order = 12
            else:
                r, a2 = self._bolza_condition_23()
                if r:
                    self.__bolza_a2 = a2
                    self.__aut_gp_order = 8
                else:
                    self.__aut_gp_order = 4

        return self.__aut_gp_order

    def mestre_U(self):
        r"""
        Return an invariant U that has a weight 12 and is nonzero in the
        context of Mestre's algorithm.

        EXAMPLES::

            sage: from pyhdme import *
            sage: PPASInvariants([2, 0, 0, 3], "Clebsch").mestre_U()
            64
            sage: PPASInvariants([0, 2, 0, 3], "Clebsch").mestre_U()
            8
            sage: PPASInvariants([0, 0, 2, 3], "Clebsch").mestre_U()
            4

        """

        if self.__mestre_U is None:
            if not self.base_ring().is_exact():
                raise ValueError("Choosing U in Mestre's algorithm is not implemented over inexact base fields")
            if self.is_geometrically_split():
                raise ValueError("Mestre's algorithm is not available for geometrically split surfaces")
            A, B, C, D = self.clebsch_invariants()
            I10 = self.igusa_clebsch_invariants()[3]
            if A != 0:
                self.__mestre_U = A**6
            elif B != 0:
                self.__mestre_U = B**3
            elif C != 0:
                self.__mestre_U = C**2
            else:
                raise ValueError("Could not find nonzero invariant of weight 12")
        return self.__mestre_U

    @cached_method
    def mestre_conic_coefficients(self):
        r"""
        Return Mestre's conic attached to the specified set of invariants as a
        set of 6 coefficients t11, t22, t33, t23, t31, t12 encoding the conic
        equation t11 * x^2 + 2 * t12 * x * y + ... + t33 * z^2 = 0. Also return
        the absolute invariants x, y, z used to compute those
        coefficients. This code is taken from
        sage.schemes.hyperelliptic_curves.mestre.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 2, 4, 4, 1]).mestre_conic_coefficients()
            ([39358/820125,
              -393120856/44840334375,
              -1213983787808/2451645281953125,
              14194881712/18160335421875,
              -393120856/44840334375,
              12083852/996451875],
             [2134/54675, 3674/2460375, -196560428/44840334375])

        """

        I2, I4, I6, I10 = self.igusa_clebsch_invariants()
        # Setting x,y,z as in Mestre's algorithm (Using Lauter and Yang's formulas)
        x = 8*(1 + 20*I4/(I2**2))/225
        y = 16*(1 + 80*I4/(I2**2) - 600*I6/(I2**3))/3375
        z = -64*(-10800000*I10/(I2**5) - 9 - 700*I4/(I2**2) + 3600*I6/(I2**3) +
                 12400*I4**2/(I2**4) - 48000*I4*I6/(I2**5))/253125
        coeffs = Sequence([x+6*y, 2*z, 6*x**2*y + 2*y**2 + 3*x*z, 9*x**3 + 4*x*y + 6*y**2, 2*z, 6*x**2+2*y],
                          universe = self.base_ring())
        # do not perform any reduction whatsoever.
        # try:
        #     den = lcm([u.denominator() for u in coeffs])
        #     coeffs = Sequence([u * den for u in coeffs], universe = self.base_ring())
        # except (AttributeError, TypeError):
        #     pass
        return coeffs, Sequence([x, y, z], universe = self.base_ring())

        # U = self.mestre_U()
        # A, B, C, D = self.clebsch_invariants()
        # I10 = self.igusa_clebsch_invariants()[3]
        # c11 = 2 * C + A * B/3
        # c22 = D
        # c33 = B * D/2 + 2 * C * (B**2 + A * C)/9
        # c23 = B * (B**2 + A * C)/3 + C * (2 * C + A * B/3)/3
        # c31 = D
        # c12 = 2 * (B**2 + A * C)/3

        # t11 = U**2 * I10**8 * c11
        # t22 = I10**10 * c22
        # t33 = U**8 * c33
        # t23 = U**4 * I10**5 * c23
        # t31 = U**5 * I10**4 * c31
        # t12 = U * I10**9 * c12
        # return Sequence([t11, t22, t33, t23, t31, t12], universe = self.base_ring())

    @cached_method
    def mestre_conic(self):
        r"""
        Return Mestre's conic attached to the specified set of invariants as a
        Conic object.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 2, 4, 4, 1]).mestre_conic().has_rational_point()
            True

        """
        c11, c22, c33, c23, c31, c12 = self.mestre_conic_coefficients()[0]
        return Conic([c11, 2 * c12, 2 * c31, c22, 2 * c23, c33])

    def mestre_line(self):
        r"""
        Return integers a2, a3, b2, b3 such that the parametric line x = t,
        y = a2 * t + b2, z = a3 * t + b3 intersects the Mestre conic in two
        distinct (non necessarily rational) points.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 2, 4, 4, 1]).mestre_line()
            [0, 0, 1, 1]

        """
        if not self.__mestre_line is None:
            return self.__mestre_line

        if not self.base_ring().is_exact():
            raise ValueError("Mestre line computation not implemented over inexact base rings")
        R = PolynomialRing(self.base_ring(), "t")
        t = R.gen()
        for i in range(5**4):
            a2 = i % 5
            a3 = (i // 5) % 5
            b2 = 1 + (i // 25) % 5
            b3 = 1 + (i // 125) % 5

            x = t
            y = a2 * t + b2
            z = a3 * t + b3
            c11, c22, c33, c23, c31, c12 = self.mestre_conic_coefficients()[0]
            substitution = c11 * x**2 + c22 * y**2 + c33 * z**2 + 2 * c23 * y * z + 2 * c31 * x * z + 2 * c12 * x * y
            c0, c1, c2 = [substitution.coefficient(i) for i in range(3)]
            delta = c1**2 - 4 * c1 * c2
            if delta != 0 and c2 != 0 and c0 != 0:
                self.__mestre_line = Sequence([a2, a3, b2, b3], universe = self.base_ring())
                break
        return self.__mestre_line

    @cached_method
    def mestre_conic_point(self):
        r"""
        Return a point over Mestre's conic defined over the base ring of self,
        raising a :class:`ValueError` if no such point can be found.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 0, 0, 2, 4, 4, 1]).mestre_conic_point()
            (102878878/211574025 : -247277/1567215 : 1)
            sage: X = PPASInvariants([1, 2, 3, 7]).change_ring(ComplexBallField(200), force_exact_computations = True)
            sage: X.mestre_conic_point()[2] == 1
            True

        """
        # do not perform any reduction whatsoever.
        try:
            pt = self.mestre_conic().rational_point()
        #     pt.clear_denominators()
        except NotImplementedError:
            a2, a3, b2, b3 = self.mestre_line()
            R = PolynomialRing(self.base_ring(), "t")
            t = R.gen()
            x = t
            y = a2 * t + b2
            z = a3 * t + b3
            c11, c22, c33, c23, c31, c12 = self.mestre_conic_coefficients()[0]
            substitution = c11 * x**2 + c22 * y**2 + c33 * z**2 + 2 * c23 * y * z + 2 * c31 * x * z + 2 * c12 * x * y
            c0, c1, c2 = [substitution.coefficient(i) for i in range(3)]
            delta = c1**2 - 4 * c0 * c2
            try:
                s = self.base_ring()(sqrt(delta))
            except (ValueError, TypeError):
                raise ValueError("Could not extract a square root of {} in {}".format(delta, base_ring))
            x = (- c1 + s) / (2 * c2)
            y = a2 * x + b2
            z = a3 * x + b3
            pt = [x, y, z]
            # try:
            #     pt = pt * pt.denominator()  # clear the denominator
            # except (AttributeError, TypeError):
            #     pass
        return pt

    def genus_2_curve_equation(self, same_invariants = False):
        r"""
        Return a genus 2 curve equation over the base ring of self which
        realizes the specified absolute invariants, raising a
        :class:`ValueError` if no such curve exists. Over inexact rings, a
        :class:`ValueError` is raised unless the invariants are the base change
        of invariants over an exact rings where some precomputations (e.g. the
        automorphism group) have been performed. The curve equation is not
        minimized.

        If "same_invariants" is False (default), the modular invariants of self
        are modified to match the computed curve equation, while keeping the
        same point in weighted projective space P(2,3,5,6). Otherwise, the
        curve equation is rescaled to match the computed invariants; this may
        fail over non-algebraically closed fields.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([1, 2, 3, 4]).absolute_igusa_invariants()
            sage: crv = PPASInvariants(vec).genus_2_curve_equation()
            sage: PPASInvariants(crv).absolute_igusa_invariants() == vec
            True
            sage: PPASInvariants(crv).modular_invariants() == [1, 2, 3, 4]
            False
            sage: vec = PPASInvariants([1, 0, 0, 2, 4, 4, 1]).modular_invariants()
            sage: crv = PPASInvariants(vec).genus_2_curve_equation(same_invariants = True)
            sage: PPASInvariants(crv).modular_invariants() == vec
            True
            sage: vec = [1, 2, 3, 4]
            sage: crv = PPASInvariants(vec).genus_2_curve_equation(same_invariants = True)
            Traceback (most recent call last):
            ...
            ValueError: Could not extract a 2th root of -776530813792251729016177857675101017435797423118579422943574066637613608874904078971127288329976286075916901/2068001658293095566095155200000000000000000000000000 in Rational Field

        TESTS::

            sage: from pyhdme import PPASInvariants
            sage: def check(coeffs): vec = PPASInvariants(coeffs).modular_invariants(); crv = PPASInvariants(vec).genus_2_curve_equation(same_invariants = True); return PPASInvariants(crv).modular_invariants() == vec
            sage: check([0, 4, 0, 0, 0, 0, 1])
            True
            sage: check([2, 0, 0, 0, 0, 0, 1])
            True
            sage: check([0, 1, 0, 0, 0, -1, 0])
            True
            sage: def check2(coeffs): vec = PPASInvariants(coeffs).absolute_igusa_invariants(); crv = PPASInvariants(vec).genus_2_curve_equation(); return PPASInvariants(crv).absolute_igusa_invariants() == vec
            sage: check2([0, 8, -12, 4, 4, -4, 1])
            True
            sage: check2([1, 4, 6, 2, 1, 2, 1])
            True

        """

        if not self.__g2_curve is None:
            return self.__g2_curve

        if self.is_geometrically_split():
            raise ValueError("The given PPAS is geometrically split")

        R = PolynomialRing(self.base_ring(), "x")
        t = R.gen()
        n = self.geometric_automorphism_group_order()

        if n == 2:
            # Mestre's algorithm
            x, y, z = self.mestre_conic_coefficients()[1]
            # setting the cijk from Mestre's algorithm
            c111 = 12*x*y - 2*y/3 - 4*z
            c112 = -18*x**3 - 12*x*y - 36*y**2 - 2*z
            c113 = -9*x**3 - 36*x**2*y - 4*x*y - 6*x*z - 18*y**2
            c122 = c113
            c123 = -54*x**4 - 36*x**2*y - 36*x*y**2 - 6*x*z - 4*y**2 - 24*y*z
            c133 = -27*x**4/2 - 72*x**3*y - 6*x**2*y - 9*x**2*z - 39*x*y**2 - \
                36*y**3 - 2*y*z
            c222 = -27*x**4 - 18*x**2*y - 6*x*y**2 - 8*y**2/3 + 2*y*z
            c223 = 9*x**3*y - 27*x**2*z + 6*x*y**2 + 18*y**3 - 8*y*z
            c233 = -81*x**5/2 - 27*x**3*y - 9*x**2*y**2 - 4*x*y**2 + 3*x*y*z - 6*z**2
            c333 = 27*x**4*y/2 - 27*x**3*z/2 + 9*x**2*y**2 + 3*x*y**3 - 6*x*y*z + \
                4*y**3/3 - 10*y**2*z
            # writing out the hyperelliptic curve polynomial
            F1, F2, F3 = PPASInvariants.parametrize_conic(self.mestre_conic_point(),
                                                          self.mestre_conic_coefficients()[0], t)
            crv = c111*F1**3 + c112*F1**2*F2 + c113*F1**2*F3 + c122*F1*F2**2 + \
                  c123*F1*F2*F3 + c133*F1*F3**2 + c222*F2**3 + c223*F2**2*F3 + \
                  c233*F2*F3**2 + c333*F3**3

            # x, y, z = PPASInvariants.parametrize_conic(self.mestre_conic_point(),
            #                                            self.mestre_conic_coefficients(), t)
            # U = self.mestre_U()
            # I10 = self.igusa_clebsch_invariants()[3]
            # A, B, C, D = self.clebsch_invariants()
            # c111 = 8 * (A**2 * C - 6 * B * C + 9 * D)/36
            # c112 = 4 * (2 * B**3 + 4 * A * B * C + 12 * C**2 + 3 * A * D)/36
            # c113 = 4 * (A * B**3 + 4 * A**2 * B * C/3 + 4 * B**2 * C + 6 * A * C**2 + 3 * B * D)/36
            # c122 = 4 * (A * B**3 + 4 * A**2 * B * C/3 + 4 * B**2 * C + 6 * A * C**2 + 3 * B * D)/36
            # c123 = 2 * (2 * B**4 + 4 * A * B**2 * C + 4 * A**2 * C**2/3 + 4 * B * C**2 + 3 * A * B * D + 12 * C * D)/36
            # c133 = 2 * (A * B**4 + 4 * A**2 * B**2 * C/3 + 16 * B**3 * C/3 + 26 * A * B * C**2/3 +  8 * C**3 + 3 * B**2 * D + 2 * A * C * D)/36
            # c222 = 4 * (3 * B**4 + 6 * A * B**2 * C + 8 * A**2 * C**2/3 + 2 * B * C**2 - 3 * C * D)/36
            # c223 = 2 * (-2 * B**3 * C/3 - 4 * A * B * C**2/3 - 4 * C**3 + 9 * B**2 * D + 8 * A * C * D)/36
            # c233 = 2 * (B**5 + 2 * A * B**3 * C + 8 * A**2 * B * C**2/9 + 2 * B**2 * C**2/3  - B * C * D + 9 * D**2)/36
            # c333 = 1 * (-2 * B**4 * C - 4 * A * B**2 * C**2 - 16 * A**2 * C**3/9 - 4 * B * C**3/3  + 9 * B**3 * D + 12 * A * B * C * D + 20 * C**2 * D)/36

            # t111 = c111 * U**3 * I10**12 * x**3
            # t112 = 3 * c112 * U**2 * I10**13 * x**2 * y
            # t113 = 3 * c113 * U**6 * I10**8 * x**2 * z
            # t122 = 3 * c122 * U * I10**14 * x * y**2
            # t123 = 6 * c123 * U**5 * I10**9 * x * y * z
            # t133 = 3 * c133 * U**9 * I10**4 * x * z**2
            # t222 = c222 * I10**15 * y**3
            # t223 = 3 * c223 * U**4 * I10**10 * y**2 * z
            # t233 = 3 * c233 * U**8 * I10**5 * y * z**2
            # t333 = c333 * U**12 * z**3
            # crv =  t111 + t112 + t113 + t122 + t123 + t133 + t222 + t223 + t233 + t333

            if same_invariants:
                # Fail if cannot find rescaling by [4, 6, 10, 12], as the only twists are quadratic twists
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                      self.minimal_weight_combination(), [4, 6, 10, 12])
                crv = crv / alpha

        elif n == 4:
            # Cardona's algorithm
            # Coefficients of conic
            A11 = 1/3*A*B + 2*C
            A12 = 2/3*B^2 + 2/3*A*C
            A22 = D
            A33 = -2/9*B^4 - 4/9*A*B^2*C - 2/9*A^2*C^2 + 1/6*A*B*D + C*D
            # Conic parametrization
            x = -2 * A22 * t - 2 * A12
            y = -A22 * t**2 + A11
            z = A22 * t**2 + 2 * A12 * t + A11
            # Substitute in cubic

            a111 = 4/675*A^2*C - 8/225*B*C + 4/75*D
            a112 = 4/675*B^3 + 8/675*A*B*C + 8/225*C^2 + 2/225*A*D
            a122 = 2/675*A*B^3 + 8/2025*A^2*B*C + 8/675*B^2*C + 4/225*A*C^2 + 2/225*B*D
            a133 = -1/2025*A^2*B^4 - 4/6075*A^3*B^2*C + 8/2025*B^5 + 14/2025*A*B^3*C + 2/2025*A^2*B*C^2 + 8/675*B^2*C^2 + 4/675*A*C^3 + 1/225*A*B^2*D + 2/675*A^2*C*D + 2/225*B*C*D - 2/75*D^2
            a222 = 2/225*B^4 + 4/225*A*B^2*C + 16/2025*A^2*C^2 + 4/675*B*C^2 - 2/225*C*D
            a233 = 1/2025*A*B^5 + 2/1215*A^2*B^3*C + 8/6075*A^3*B*C^2 - 2/2025*B^4*C + 2/2025*A*B^2*C^2 + 8/2025*A^2*C^3 - 4/675*B*C^3 + 2/675*B^3*D + 1/675*A*B*C*D - 2/225*C^2*D - 1/225*A*D^2

            t111 = -A33 * a111 * x**3
            t112 = -3 * A33 * a112 * x**2 * y
            t122 = -3 * A33 * a122 * x * y**2
            t133 = 3 * A22 * a133 * x * z**2
            t222 = -A33 * a222 * y**3
            t333 = 3 * A22 * a233 * x * z**2
            crv = t111 + t112 + t122 + t133 + t222 + t333
            if same_invariants:
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                try:
                    alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                          self.minimal_weight_combination(), [4, 6, 10, 12])
                    crv = crv / alpha
                except ValueError:
                    raise NotImplementedError("Curve reconstruction with the same invariants over a non-algebraically closed field in the presence of geometric automorphisms is not fully implemented")

        elif n == 8:
            crv =  t**5 + t**3 + (1 / self._bolza_a2()) * t
            if same_invariants:
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                try:
                    alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                          self.minimal_weight_combination(), [4, 6, 10, 12])
                    crv = crv / alpha
                except ValueError:
                    raise NotImplementedError("Curve reconstruction with the same invariants over a non-algebraically closed field in the presence of geometric automorphisms is not fully implemented")

        elif n == 10:
            crv = t**6 + t
            if same_invariants:
                new_IC = crv.discriminant()
                alpha = PPASInvariants.find_rescaling(self.base_ring(), [new_IC], [self.igusa_clebsch_invariants()[3]],
                                                      [1], [2])
                crv = alpha * t**6 + t / alpha

        elif n == 12:
            crv =  t**6 + t**3 + (1 / self._bolza_a2())
            if same_invariants:
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                try:
                    alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                          self.minimal_weight_combination(), [4, 6, 10, 12])
                    crv = crv / alpha
                except ValueError:
                    raise NotImplementedError("Curve reconstruction with the same invariants over a non-algebraically closed field in the presence of geometric automorphisms is not fully implemented")

        elif n == 24:
            crv = t**6 + 1
            if same_invariants:
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                      self.minimal_weight_combination(), [2, 3, 5, 6])
                crv = t**6 + (1 / alpha)

        else:
            assert n == 48, "Unknown automorphism group order: {}".format(n)
            crv = t**5 + t
            if same_invariants:
                new_IC = PPASInvariants.modified_igusa_clebsch_from_curve(crv)
                alpha = PPASInvariants.find_rescaling(self.base_ring(), new_IC, self.modified_igusa_clebsch_invariants(),
                                                      self.minimal_weight_combination(), [2, 3, 5, 6])
                crv = t**5 + (1 / alpha) * t

        # Make coefficients look nicer?
        # try:
        #     crv = crv * crv.denominator()
        #     crv = crv / gcd([x.numerator() for x in crv.coefficients()])
        # except (AttributeError, TypeError):
        #     pass
        # Adjust genus 2 curve equation to achieve the given invariants
        # Make coefficients look nicer
        # if self.base_ring().is_exact():
        #     u = crv.coefficient(4)
        #     v = crv.coefficient(2)
        #     if u != 0:
        #         crv = u**3 * crv.subs(t / u)
        #     elif v != 0:
        #         crv = v**(-3) * crv.subs(t * v)

        self.__g2_curve = crv
        return self.__g2_curve

    def j_invariants_quadratic_equation(self):
        r"""
        Return a quadratic equation over the base field of self whose solutions
        are the j-invariants of the two elliptic factors of the specified PPAS.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: eq = PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).j_invariants_quadratic_equation()
            sage: t = eq.parent().gen()
            sage: eq == (t - 432/7) * (t - 1728/5)
            True

        """
        if not self.is_geometrically_split():
            raise ValueError("The given PPAS is not geometrically split")

        X4, X6, X10, X12 = self.modular_invariants()
        j1j2 = 12 * X4**3 / X12
        j1pj2 = (12 / X12) * ((X4**3 - X6**2 + 12**5 * X12) / 1728)
        R = PolynomialRing(self.base_ring(), "t")
        t = R.gen()
        return t**2 - j1pj2 * t + j1j2

    def elliptic_curves(self):
        r"""
        Return the two elliptic factors of the PPAS with the specified
        invariants. We require that they are defined over the base field.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: vec = PPASInvariants([EllipticCurve([1, 2]), EllipticCurve([3, 4])]).modular_invariants()
            sage: E1, E2 = PPASInvariants(vec).elliptic_curves()
            sage: {E1.j_invariant(), E2.j_invariant()} == {432/7, 1728/5}
            True
            sage: PPASInvariants([E1, E2]).modular_invariants() == vec
            True

        """

        if self.__ell_curves is None:
            L = self.j_invariants_quadratic_equation().roots(multiplicities = False)
            if len(L) == 0:
                raise ValueError("Elliptic factors are not defined over the base field")
            elif len(L) == 1:
                E2 = EllipticCurve_from_j(L[0])
            else:
                E2 = EllipticCurve_from_j(L[1])
            E1 = EllipticCurve_from_j(L[0])
            c41, c61 = E1.c_invariants()
            c42, c62 = E2.c_invariants()
            delta1 = E1.discriminant()
            delta2 = E2.discriminant()
            new_invs = [c41 * c42, c61 * c62, 0, (delta1 * delta2) / 12]
            u = PPASInvariants.find_rescaling(self.base_ring(), new_invs, self.modular_invariants(),
                                              self.minimal_weight_combination(), [2, 3, 5, 6])
            _ , _, _, a4, a6 = E1.a_invariants()
            E1 = EllipticCurve([0, 0, 0, a4 / u**2, a6 / u**3])
            self.__ell_curves = [E1, E2]

        return self.__ell_curves

    def change_ring(self, R, force_exact_computations = False):
        r"""
        Return an object of the PPASInvariants class with the same invariants
        over the ring R. If `force_exact_computations` is set to True, then all
        computations that cannot be performed over an inexact ring are
        performed before the base change.

        EXAMPLES::

            sage: from pyhdme import PPASInvariants
            sage: PPASInvariants([1, 2, 13, 4]).change_ring(FiniteField(13)).is_geometrically_split()
            True
            sage: X = PPASInvariants([1, 0, 0, 0, 0, 0, 1/3]).change_ring(ComplexBallField(500), force_exact_computations = True)
            sage: X.geometric_automorphism_group_order()
            24
            sage: X.is_geometrically_split()
            False

        """

        R = PPASInvariants.ambient_field(R)
        reduction = (self.base_ring().characteristic() == 0) and (R.characteristic() > 0)

        if force_exact_computations and not reduction:
            if not self.is_geometrically_split():
                self.geometric_automorphism_group_order()
                self.mestre_U()
                self.mestre_line()
                self.minimal_weight_combination()
            else:
                pass

        res = PPASInvariants(Sequence(self.modular_invariants(), universe = R))
        if not self.__g2_curve is None and not reduction:
            res.__g2_curve = self.__g2_curve.change_ring(R)
        if not self.__ell_curves is None:
            res.__ell_curves = [E.base_extend(R) for E in self.__ell_curves]
        if not self.__ic is None:
            res.__ic = Sequence(self.__ic, universe = R)
        if not self.__clebsch is None:
            res.__clebsch = Sequence(self.__clebsch, universe = R)
        if not self.__ic_mod is None:
            res.__ic_mod = Sequence(self.__ic_mod, universe = R)
        if not self.__mestre_U is None and not reduction:
            res.__mestre_U = R(self.__mestre_U)
        if not self.__bolza_a2 is None and not reduction:
            res.__bolza_a2 = R(self.__bolza_a2)
        if not self.__mestre_line is None and not reduction:
            res.__mestre_line = list(self.__mestre_line)
        if not self.__aut_gp_order is None and not reduction:
            res.__aut_gp_order = self.__aut_gp_order
        if not self.__min_wt is None and not reduction:
            res.__min_wt = list(self.__min_wt)
        return res
