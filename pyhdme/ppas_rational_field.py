
class PPAS_rational_field(SageObject):

    def __init__(self, data, inv_type = "IgusaClebsch"):

        self.invariants = PPASInvariants.__init__(data, inv_type)

        self.cplx_prec = None
        self.cplx_g2_curve = None
        self.cplx_elliptic_curves = None
        self.curve_discrete_data = None
        self.periods = None
        self.cplx_base_change = None

