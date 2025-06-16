
class RigidifiedPPAS(SageObject):

    def _check_base_ring(F):
        if not F.is_ring():
            raise TypeError("Base ring must be a ring")
        try:
            F(1)/2
            F(1)/3
            F(1)/5
        except:
            raise ValueError("2, 3 and 5 must be invertible in the base ring")

    def __init__(self, data, inv_type = "IgusaClebsch"):

        self.base_ring = None
        self.g2_curve = None
        self.elliptic_curves = None

        self.mfgens = None
        self.ic = None
        self.clebsch = None
        self.ic_mod = None
        self.igusa_integral = None

        self.cplx_prec = None
        self.cplx_curve = None
        self.cplx_base_change = None
        self.cplx_base_change_det = None
        self.periods = None

        if isinstance(data, list):

            if len(data) == 2:
                # Interpret data as a pair of elliptic curves
                E1 = data[0]
                E2 = data[1]
                F = Sequence(E1.a_invariants() + E2.a_invariants()).universe()
                RigidifiedPPAS._check_base_ring(F)
                self.base_ring = F
                self.elliptic_curves = [E1.change_ring(F), E2.change_ring(F)]
                c41, c61 = E1.c_invariants()
                c42, c62 = E2.c_invariants()
                delta1 = E1.discriminant()
                delta2 = E2.discriminant()
                self.mfgens = Sequence([c41 * c42, c61 * c62, 0, delta1 * delta2], universe = F)

            elif len(data) == 3:
                # Interpret data as a triple of three absolute Igusa invariants I4*I6'/I10, I4^2*I12/I10^2, I4^5/I10^2
                F = Sequence(data).universe()
                RigidifiedPPAS._check_base_ring(F)
                self.base_ring = F

                self.ic_mod = Sequence([data[2], data[0] * data[2], data[2]**2, data[1] * data[2]**2], universe = F)
                self.ic = g2_curve_ic_from_ic_mod(self.ic_mod)
                self.mfgens = g2_curve_mfgens_from_ic(self.ic)

            elif len(data) == 4:
                # Interpret data as a tuple of invariants of the specified type
                F = Sequence(data).universe()
                RigidifiedPPAS._check_base_ring(F)
                self.base_ring = F
                data = Sequence(data, universe = F)

                if inv_type == "IgusaClebsch":
                    self.ic = data
                    self.mfgens = g2_curve_mfgens_from_ic(self.ic)
                elif inv_type == "Clebsch":
                    self.clebsch = data
                    self.ic = g2_curve_ic_from_clebsch(self.clebsch)
                    self.mfgens = g2_curve_mfgens_from_ic(self.ic)
                elif inv_type == "ModifiedIgusaClebsch":
                    self.ic_mod = data
                    self.ic = g2_curve_ic_from_ic_mod(self.ic_mod)
                    self.mfgens = g2_curve_mfgens_from_ic(self.ic)
                elif inv_type == "ModularForms":
                    self.mfgens = data
                else:
                    raise TypeError("Unknown invariant type: {}".format(inv_type))

            elif len(data) == 7:
                # Interpret data as the coefficients of a genus 2 curve
                F = Sequence(data).universe()
                RigidifiedPPAS._check_base_ring(F)
                self.base_ring = F

                R = PolynomialRing(F, "x")
                self.g2_curve = R(data).reverse(6)
                self.mfgens = g2_curve_mfgens(self.g2_curve)

            else:
                raise TypeError("Invalid input length {}".format(len(data)))

        else:
            # Interpret data as an element of a polynomial ring
            pass
