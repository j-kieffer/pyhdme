r"""
Invariants and reconstruction of principally polarized abelian surfaces

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

# Copyright 2025 Jean Kieffer
# See LICENSE file for license details.

class PPASInvariants(SageObject):
    r"""
    Create a data structure allowing conversions between different types of
    invariants of principally polarized abelian surfaces over a field.

    We require that 2, 3 and 5 are invertible in the base field.

    The different sets of invariants are as follows:

    1. The "modular invariants", which correspond to the following Siegel
    modular forms with integral Fourier expansions:

    psi4 = 1 + 240(q1+q2) + ...
    psi6 = 1 - 504(q1+q1) + ...
    chi10 = (q3 - 2 + q3^-1) + ...
    chi12 = (q3 + 10 + q3^-1) + ...

    These invariants make sense both for Jacobians of genus 2 curves and for
    products of elliptic curves.

    2. The classical Igusa--Clebsch invariants I_2, I_4, I_6, I_{10}.

    3. The classical Clebsch invariants A, B, C, D.

    4. The modified Igusa--Clebsch invariants I_4, I_6', I_{10}, I_{12}.

    5. The absolute Igusa invariants, defined as follows:
    j_1 = I_4*I_6/I_{10}, j_2 = I_4^2*I_{12}/I_{10}^2, j_3 = I_4^5/I_{10}^2.

    6. The equation of a genus 2 curve y^2 = f(x), represented as the
    polynomial f(x) of degree 5 or 6.

    7. A pair of j-invariants of elliptic curves.

    8. A pair of equations of elliptic curves.

    An element of the PPASInvariants class may be initialized using any of
    these types of invariants, and supports various conversions. A
    :class:`ValueError` is raised when the conversion doesn't make sense
    (e.g. when asking for the Igusa invariants of a product of elliptic curves)

    Some properties of principally polarized abelian surfaces that are easily
    obtained from invariants (automorphism groups) are also available.

    EXAMPLES::

        sage:

    TESTS::

        sage:

    """

    def field_base_change(F):
        r"""
        Return a minimal base field F' that contains the ring F, and raises
        a :class:`ValueError` if F' has characteristic 2, 3, or 5.

        EXAMPLES::

            sage:

        TESTS::

            sage:

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

            sage:

        TESTS::

            sage:

        """
        I4, I6p, I10, I12 = vec
        return [I4 / 4, I6p / 4, -I10 / 2**12, I12 / 2**15]

    def modified_igusa_clebsch_from_igusa_clebsch(vec):
        r"""
        Return the modified Igusa--Clebsch invariants associated to the given Igusa--Clebsch invariants.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        I2, I4, I6, I10 = vec
        return [I4, (I2 * I4 - 3 * I6) / 2, I10, I2 * I10]

    def igusa_clebsch_from_clebsch(vec):
        r"""
        Return the Igusa--Clebsch invariants associated to the given Clebsch invariants.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        A, B, C, D = vec
        I2 = - 120 * A
        I4 = - 720 * A**2 + 6750 * B
        I6 = 8640 * A**3 - 108000 * A * B + 202500 * C
        I10 = - 62208 * A**5 + 972000 * A**3 * B + 1620000 * A**2 * C - 3037500 * A * B**2 - 6075000 * B * C - 4556250 * D
        return [I2, I4, I6, I10]

    def modified_igusa_clebsch_from_curve_coefficients(vec):
        r"""
        Return the modified Igusa--Clebsch invariants from the given vector
        of curve coefficients. Todo: rewrite this in terms of transvectants
        of binary forms for efficiency.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        a, b, c, d, e, f, g = vec
        I2 = -240 * g * a + (40 * f * b + (-16 * e * c + 6 * d**2))
        I4 = 1620 * g**2 * a**2 + (-540 * g * f * b + ((-504 * g * e + 300 * f**2) * c + (324 * g * d**2 - 180 * f * e * d + 48 * e**3))) * a + ((300 * g * e - 80 * f**2) * b**2 + ((-180 * g * d + 4 * f * e) * c + (36 * f * d**2 - 12 * e**2 * d)) * b + (48 * g * c**3 + (-12 * f * d + 4 * e**2) * c**2))
        I6p = -14580 * g**3 * a**3 + (7290 * g**2 * f * b + ((16524 * g**2 * e -8100 * g * f**2) * c + (-18954 * g**2 * d**2 + (17010 * g * f * e - 3375 * f**3) * d+ (-5616 * g * e**3 + 1350 * f**2 * e**2)))) * a**2 + ((-8100 * g**2 * e +2160 * g * f**2) * b**2 + ((17010 * g**2 * d + (-11448 * g * f * e +3600 * f**3)) * c + (-2187 * g * f * d**2 + (2754 * g * e**2 - 810 * f**2 * e) * d +36 * f * e**3)) * b + (-5616 * g**2 * c**3 + (2754 * g * f * d + (2916 * g * e**2 -1440 * f**2 * e)) * c**2 + ((-3402 * g * e + 405 * f**2) * d**2 + 702 * f * e**2 * d- 144 * e**4) * c + (729 * g * d**4 - 243 * f * e * d**3 + 54 * e**3 * d**2))) * a +((-3375 * g**2 * d + (3600 * g * f * e - 1120 * f**3)) * b**3 + (1350 * g**2 * c**2+ (-810 * g * f * d + (-1440 * g * e**2 + 624 * f**2 * e)) * c + ((405 * g * e +216 * f**2) * d**2 - 279 * f * e**2 * d + 54 * e**4)) * b**2 + (36 * g * f * c**3 +((702 * g * e - 279 * f**2) * d + 6 * f * e**2) * c**2 + (-243 * g * d**3 +81 * f * e * d**2 - 18 * e**3 * d) * c) * b + ((-144 * g * e + 54 * f**2) * c**4 +(54 * g * d**2 - 18 * f * e * d + 4 * e**3) * c**3))
        I10 = -46656 * g**5 * a**5 + (38880 * g**4 * f * b + ((62208 * g**4 * e -32400 * g**3 * f**2) * c + (34992 * g**4 * d**2 + (-77760 * g**3 * f * e +27000 * g**2 * f**3) * d + (-13824 * g**3 * e**3 + 43200 * g**2 * f**2 * e**2 -22500 * g * f**4 * e + 3125 * f**6)))) * a**4 + ((-32400 * g**4 * e +540 * g**3 * f**2) * b**2 + ((-77760 * g**4 * d + (31968 * g**3 * f * e -1800 * g**2 * f**3)) * c + (15552 * g**3 * f * d**2 + (46656 * g**3 * e**2 -31320 * g**2 * f**2 * e + 2250 * g * f**4) * d + (-21888 * g**2 * f * e**3 +15600 * g * f**3 * e**2 - 2500 * f**5 * e))) * b + (-13824 * g**4 * c**3 +(46656 * g**3 * f * d + (-17280 * g**3 * e**2 - 6480 * g**2 * f**2 * e +1500 * g * f**4)) * c**2 + ((3888 * g**3 * e - 27540 * g**2 * f**2) * d**2 +(-3456 * g**2 * f * e**2 + 19800 * g * f**3 * e - 3750 * f**5) * d +(9216 * g**2 * e**4 - 10560 * g * f**2 * e**3 + 2000 * f**4 * e**2)) * c +(-8748 * g**3 * d**4 + (21384 * g**2 * f * e - 1350 * g * f**3) * d**3 +(-8640 * g**2 * e**3 - 9720 * g * f**2 * e**2 + 2250 * f**4 * e) * d**2 +(6912 * g * f * e**4 - 1600 * f**3 * e**3) * d + (-1024 * g * e**6 +256 * f**2 * e**5)))) * a**3 + ((27000 * g**4 * d + (-1800 * g**3 * f * e +410 * g**2 * f**3)) * b**3 + (43200 * g**4 * c**2 + (-31320 * g**3 * f * d +(-6480 * g**3 * e**2 + 8748 * g**2 * f**2 * e - 1700 * g * f**4)) * c +((-27540 * g**3 * e + 15417 * g**2 * f**2) * d**2 + (16632 * g**2 * f * e**2 -12330 * g * f**3 * e + 2000 * f**5) * d + (-192 * g**2 * e**4 + 248 * g * f**2 * e**3- 50 * f**4 * e**2))) * b**2 + (-21888 * g**3 * f * c**3 + ((-3456 * g**3 * e +16632 * g**2 * f**2) * d + (15264 * g**2 * f * e**2 - 13040 * g * f**3 * e +2250 * f**5)) * c**2 + (21384 * g**3 * d**3 + (-22896 * g**2 * f * e +1980 * g * f**3) * d**2 + (-5760 * g**2 * e**3 + 10152 * g * f**2 * e**2 -2050 * f**4 * e) * d + (-640 * g * f * e**4 + 160 * f**3 * e**3)) * c +(-6318 * g**2 * f * d**4 + (5832 * g**2 * e**2 + 3942 * g * f**2 * e -900 * f**4) * d**3 + (-4464 * g * f * e**3 + 1020 * f**3 * e**2) * d**2 +(768 * g * e**5 - 192 * f**2 * e**4) * d)) * b + ((9216 * g**3 * e -192 * g**2 * f**2) * c**4 + (-8640 * g**3 * d**2 + (-5760 * g**2 * f * e -120 * g * f**3) * d + (-4352 * g**2 * e**3 + 4816 * g * f**2 * e**2 -900 * f**4 * e)) * c**3 + (5832 * g**2 * f * d**3 + (8208 * g**2 * e**2 -4536 * g * f**2 * e + 825 * f**4) * d**2 + (-2496 * g * f * e**3 +560 * f**3 * e**2) * d + (512 * g * e**5 - 128 * f**2 * e**4)) * c**2 +((-4860 * g**2 * e + 162 * g * f**2) * d**4 + (2808 * g * f * e**2 -630 * f**3 * e) * d**3 + (-576 * g * e**4 + 144 * f**2 * e**3) * d**2) * c +(729 * g**2 * d**6 + (-486 * g * f * e + 108 * f**3) * d**5 + (108 * g * e**3 -27 * f**2 * e**2) * d**4))) * a**2 + ((-22500 * g**4 * c + (2250 * g**3 * f * d +(1500 * g**3 * e**2 - 1700 * g**2 * f**2 * e + 320 * g * f**4))) * b**4 +(15600 * g**3 * f * c**2 + ((19800 * g**3 * e - 12330 * g**2 * f**2) * d +(-13040 * g**2 * f * e**2 + 9768 * g * f**3 * e - 1600 * f**5)) * c +(-1350 * g**3 * d**3 + (1980 * g**2 * f * e - 208 * g * f**3) * d**2 +(-120 * g**2 * e**3 - 682 * g * f**2 * e**2 + 160 * f**4 * e) * d + (144 * g * f * e**4- 36 * f**3 * e**3))) * b**3 + ((-10560 * g**3 * e + 248 * g**2 * f**2) * c**3 +(-9720 * g**3 * d**2 + (10152 * g**2 * f * e - 682 * g * f**3) * d +(4816 * g**2 * e**3 - 5428 * g * f**2 * e**2 + 1020 * f**4 * e)) * c**2 +(3942 * g**2 * f * d**3 + (-4536 * g**2 * e**2 - 2412 * g * f**2 * e +560 * f**4) * d**2 + (3272 * g * f * e**3 - 746 * f**3 * e**2) * d + (-576 * g * e**5+ 144 * f**2 * e**4)) * c + (162 * g**2 * e * d**4 + (-108 * g * f * e**2 +24 * f**3 * e) * d**3 + (24 * g * e**4 - 6 * f**2 * e**3) * d**2)) * b**2 +((6912 * g**3 * d + (-640 * g**2 * f * e + 144 * g * f**3)) * c**4 +(-4464 * g**2 * f * d**2 + (-2496 * g**2 * e**2 + 3272 * g * f**2 * e -630 * f**4) * d + (-96 * g * f * e**3 + 24 * f**3 * e**2)) * c**3 + ((2808 * g**2 * e- 108 * g * f**2) * d**3 + (-1584 * g * f * e**2 + 356 * f**3 * e) * d**2 +(320 * g * e**4 - 80 * f**2 * e**3) * d) * c**2 + (-486 * g**2 * d**5 +(324 * g * f * e - 72 * f**3) * d**4 + (-72 * g * e**3 +18 * f**2 * e**2) * d**3) * c) * b + (-1024 * g**3 * c**6 + (768 * g**2 * f * d +(512 * g**2 * e**2 - 576 * g * f**2 * e + 108 * f**4)) * c**5 + ((-576 * g**2 * e +24 * g * f**2) * d**2 + (320 * g * f * e**2 - 72 * f**3 * e) * d + (-64 * g * e**4 +16 * f**2 * e**3)) * c**4 + (108 * g**2 * d**4 + (-72 * g * f * e + 16 * f**3) * d**3 +(16 * g * e**3 - 4 * f**2 * e**2) * d**2) * c**3)) * a + (3125 * g**4 * b**6 +(-2500 * g**3 * f * c + ((-3750 * g**3 * e + 2000 * g**2 * f**2) * d +(2250 * g**2 * f * e**2 - 1600 * g * f**3 * e + 256 * f**5))) * b**5 +((2000 * g**3 * e - 50 * g**2 * f**2) * c**2 + (2250 * g**3 * d**2 +(-2050 * g**2 * f * e + 160 * g * f**3) * d + (-900 * g**2 * e**3 +1020 * g * f**2 * e**2 - 192 * f**4 * e)) * c + (-900 * g**2 * f * d**3 +(825 * g**2 * e**2 + 560 * g * f**2 * e - 128 * f**4) * d**2 + (-630 * g * f * e**3 +144 * f**3 * e**2) * d + (108 * g * e**5 - 27 * f**2 * e**4))) * b**4 +((-1600 * g**3 * d + (160 * g**2 * f * e - 36 * g * f**3)) * c**3 +(1020 * g**2 * f * d**2 + (560 * g**2 * e**2 - 746 * g * f**2 * e + 144 * f**4) * d +(24 * g * f * e**3 - 6 * f**3 * e**2)) * c**2 + ((-630 * g**2 * e +24 * g * f**2) * d**3 + (356 * g * f * e**2 - 80 * f**3 * e) * d**2 + (-72 * g * e**4 +18 * f**2 * e**3) * d) * c + (108 * g**2 * d**5 + (-72 * g * f * e + 16 * f**3) * d**4+ (16 * g * e**3 - 4 * f**2 * e**2) * d**3)) * b**3 + (256 * g**3 * c**5 +(-192 * g**2 * f * d + (-128 * g**2 * e**2 + 144 * g * f**2 * e - 27 * f**4)) * c**4 +((144 * g**2 * e - 6 * g * f**2) * d**2 + (-80 * g * f * e**2 + 18 * f**3 * e) * d +(16 * g * e**4 - 4 * f**2 * e**3)) * c**3 + (-27 * g**2 * d**4 + (18 * g * f * e -4 * f**3) * d**3 + (-4 * g * e**3 + f**2 * e**2) * d**2) * c**2) * b**2)
        return [I4, I6p, I10, I2 * I10]

    def __init__(self, data, inv_type = "Modular"):
        r"""
        Initialize a PPASInvariants data structure.

        The input data can be either:
        - a polynomial f of degree 5 or 6 encoding the genus 2 curve y^2 = f(x),
        - a pair of elliptic curves,
        - a triple of absolute Igusa invariants,
        - a tuple of 4 invariants, whose type is specified by inv_type. The
          possible values are: Modular (default), IgusaClebsch,
          ModifiedIgusaClebsch, and Clebsch, or
        - a tuple of 7 coefficients a_6, ..., a_0, encoding the polynomial
          f = a_6 x^6 + ... + a_0 as in the first item.

        The corresponding modular invariants are then computed.

        EXAMPLES::

            sage:

        """

        self.base_ring = None
        self.g2_curve = None
        self.elliptic_curves = None

        self.modular = None
        self.ic = None
        self.clebsch = None
        self.ic_mod = None
        self.j_invariants = None

        self.aut_gp = None

        if isinstance(data, CommutativePolynomial):
            F = data.parent().base_ring()
            F = PPASInvariants.field_base_change(F)
            self.base_ring = F
            self.g2_curve = data.base_extend(F)
            self.ic_mod = Sequence(PPASInvariants.modified_igusa_from_curve_coefficients(self.g2_curve.coefficients(6)),
                                   universe = F)
            self.modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.ic_mod),
                                   universe = F)

        elif isinstance(data, list):

            if len(data) == 2:
                E1 = data[0]
                E2 = data[1]
                if not (isinstance(E1, EllipticCurve_generic) and isinstance(E2, EllipticCurve_generic)):
                    raise TypeError("Input must be a pair of elliptic curves")
                F = Sequence(E1.a_invariants() + E2.a_invariants()).universe()
                F = PPASInvariants.field_base_change(F)
                self.base_ring = F
                self.elliptic_curves = [E1.change_ring(F), E2.change_ring(F)]
                c41, c61 = E1.c_invariants()
                c42, c62 = E2.c_invariants()
                delta1 = E1.discriminant()
                delta2 = E2.discriminant()
                self.modular = Sequence([c41 * c42, c61 * c62, 0, delta1 * delta2], universe = F)

            elif len(data) == 3:
                F = Sequence(data).universe()
                F = PPASInvariants.field_base_change(F)
                self.base_ring = F
                self.ic_mod = Sequence([data[2], data[0] * data[2], data[2]**2, data[1] * data[2]**2], universe = F)
                self.modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.ic_mod),
                                       universe = F)

            elif len(data) == 4:
                F = Sequence(data).universe()
                F = PPASInvariants.field_base_change(F)
                self.base_ring = F
                data = Sequence(data, universe = F)

                if inv_type == "Clebsch":
                    self.clebsch = data
                    data = Sequence(PPASInvariants.igusa_clebsch_from_clebsch(self.clebsch), universe = F)
                if inv_type in ["Clebsch", "IgusaClebsch"]:
                    self.ic = data
                    data = Sequence(PPASInvariants.modified_igusa_clebsch_from_igusa_clebsch(self.ic),
                                    universe = F)
                if inv_type in ["Clebsch", "IgusaClebsch", "ModifiedIgusaClebsch"]:
                    self.ic_mod = data
                    self.modular = Sequence(PPASInvariants.modular_from_modified_igusa_clebsch(self.ic_mod),
                                            universe = F)
                elif inv_type == "Modular":
                    self.modular = data
                else:
                    raise TypeError("Unknown invariant type: {}".format(inv_type))

            elif len(data) == 7:
                F = Sequence(data).universe()
                F = PPASInvariants.field_base_change(F)
                self.base_ring = F

                R = PolynomialRing(F, "x")
                self.g2_curve = R(data).reverse(6)
                self.ic_mod = Sequence(PPASInvariants.modified_igusa_from_curve_coefficients(self.g2_curve.coefficients(6)),
                                       universe = F)
                self.modular = Sequence(PPASInvariants.modular_from_modified_igusa(self.ic_mod),
                                        universe = F)

            else:
                raise TypeError("Invalid input length {}".format(len(data)))


        else:
            raise TypeError("Input must be either a polynomial, a pair of elliptic curves, or a list of coefficients or invariants")

    def base_ring(self):
        r"""
        Return the base ring of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        return self.base_ring

    def modular_invariants(self):
        r"""
        Return the modular invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        return self.modular

    def modified_igusa_clebsch_invariants(self):
        r"""
        Return the modified Igusa--Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        if self.ic_mod is None:
            m4, m6, m10, m12 = self.modular
            self.ic_mod = Sequence([4 * m4, 4 * m6, - 2**12 * m10, 2**15 * m12], universe = self.base_ring)
        return self.ic_mod

    def is_geometrically_split(self):
        r"""
        Return True iff the provided invariants correspond to a geometrically split PPAS.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        return self.modular[2] == 0

    def igusa_clebsch_invariants(self):
        r"""
        Return the Igusa--Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        if self.ic is None:
            if self.is_geometrically_split():
                raise ValueError("Igusa--Clebsch or Clebsch invariants are not defined for products of elliptic curves")
            I4, I6p, I10, I12 = self.modified_igusa_clebsch_invariants()
            I2 = I12 / I10
            I6 = (2 * I6p - I2 * I4) / (-3)
            self.ic = Sequence([I2, I4, I6, I10], universe = self.base_ring)
        return self.ic

    def clebsch_invariants(self):
        r"""
        Return the Clebsch invariants of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        if self.clebsch is none:
            I2, I4, I6, I10 = self.igusa_clebsch_invariants()
            A = -I2 / 120
            B = (I4 + 720 * I2**2) / 6750
            C = (I6 - 8640 * I2**3 + 108000 * I2 * I4) / 202500
            D = (I10 + 62208 * I2**5 - 972000 * I2**3 * I4 - 1620000 * I2**2 * I6 + 3037500 * I2 * I4**2 + 6075000 * I4 * I6) / 4556250
            self.clebsch = Sequence([A, B, C, D], universe = F)
        return self.clebsch

    @cached
    def R2_invariant(self):
        r"""
        Return the invariant R^2, of weight 30, of the principally polarized abelian surface.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        a, b, c, d = self.igusa_clebsch_invariants()
        R2 = 125971200000 * d**3 + 236196 * d**2 * a**5 + 19245600 * d**2 * a**3 * b - 507384000 * d**2 * a * b**2 - 972 * d * a**6 * b**2 - 77436 * d * a**4 * b**3 + 592272 * d * a**2 * b**4 + a**7 * b**4 - 41472 * d * b**5 + 78 * a**5 * b**5 - 159 * a**3 * b**6 + 80 * a * b**7 - 104976000 * d**2 * a**2 * c + 2099520000 * d**2 * b * c + 5832 * d * a**5 * b * c + 870912 * d * a**3 * b**2 * c - 4743360 * d * a * b**3 * c - 12 * a**6 * b**3 * c - 1332 * a**4 * b**4 * c + 1728 * a**2 * b**5 * c - 384 * b**6 * c - 8748 * d * a**4 * c**2 - 3090960 * d * a**2 * b * c**2 + 9331200 * d * b**2 * c**2 + 54 * a**5 * b**2 * c**2 + 8910 * a**3 * b**3 * c**2 - 6048 * a * b**4 * c**2 + 3499200 * d * a * c**3 - 108 * a**4 * b * c**3 - 29376 * a**2 * b**2 * c**3 + 6912 * b**3 * c**3 + 81 * a**3 * c**4 + 47952 * a * b * c**4 - 31104 * c**5
        return R2

    @cached
    def igusa_integral_invariants(self):
        r"""
        Return the values of the 14 generators X_4, X_6, ..., X_{48} for the
        ring of even-weight Siegel modular forms over `\ZZ`, of weights 4, 6,
        10, 12, 12, 16, 18, 24, 28, 30, 36, 40, 42, 48 for the given PPAS.

        EXAMPLES::

            sage:

        TESTS::

            sage:

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
                        universe = self.base_ring)

    def bolza_condition_19(self):
        r"""
        Returns True iff Bolza's condition 19 holds. If so, additionally return
        a^2 such that the given PPAS is isomorphic to the Jacobian of the genus
        2 curve y^2 = x^6 + a*x^3 + 1.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        A, B, C, D = self.clebsch_invariants()
        t1 = - B**3 + 6 * C**2
        t4 = A * B + 6 * C
        t2 = -2 * t4 * B + 9 * D
        t3 = -15 * C + 2 * A * B
        r = (t1 == 0) and (t2 == 0)
        if r:
            if (D == 0 || t3 == 0):
                raise ValueError("Unexpected vanishing")
            a2 = 100 * (-t4 / t3)
        else:
            a2 = None
        return r, a2

    def bolza_condition_23(self):
        r"""
        Returns True iff Bolza's condition 23 holds. If so, additionally return
        a^2 such that the given PPAS is isomorphic to the Jacobian of the genus
        2 curve y^2 = x^5 + a*x^3 + x.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        A, B, C, D = self.clebsch_invariants()
        t1 = B**2 * A**3 - 6 * B * C + 4 * C * A**2 - 18 * D
        t2 = 4 * B**3 + 5 * C * B * A + 6 * C**2 - 3 * A * D
        t3 = 6 * C**2 - B**3
        r = (t1 == 0) and (t2 == 0)
        if r:
            if (D == 0 || t3 == 0):
                raise ValueError("Unexpected vanishing")
            a2 = 100 * (B**2 + A * C) / (2 * A**2 * B - 3 * A * C - 15 * B**2)
        else:
            a2 = None
        return r, a2

    def abstract_automorphism_group(self):
        r"""
        Return the geometric automorphism group of the specified principally
        polarized abelian surface as an abstract group.

        EXAMPLES::

            sage:

        TESTS::

            sage:

        """
        if not self.aut_gp is None:
            return self.aut_gp

        if not self.base_ring.is_exact():
            raise NotImplementedError("Automorphism groups are not currently implemented over inexact fields")
        if self.is_geometrically_split():
            raise NotImplementedError("Automorphism groups are not currently implemented for products of elliptic curves")

        R2 = self.R2_invariant()
        A, B, C, D = self.clebsch_invariants()
        if R2 != 0 and (A != 0 or B != 0 or C != 0):
            return AbelianGroup([2])
        elif R2 != 0:
            pass
        elif B == 0 and C == 0 and D == 0:
            pass
        elif 6 * B - A**2 == 0 and 6 * C - A * B == 0 and D == 0:
            pass
        elif self.bolza_condition_19():
            pass
        elif self.bolza_condition_23():
            pass
        else:
            pass

    def genus_2_curve_equation(self):
        if self.g2_curve is None:
            if self.is_geometrically_split():
                raise ValueError("The given PPAS is geometrically split")
            self.g2_curve = g2_curve_from_igusa_clebsch(self.igusa_clebsch_invariants())
        return self.g2_curve

    def elliptic_curves(self):
        if self.elliptic_curves is None:
            if not self.is_geometrically_split():
                raise ValueError("The given PPAS is not geometrically split")
            pass
        return self.elliptic_curves
