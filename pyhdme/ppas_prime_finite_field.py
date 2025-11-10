r"""
Principally polarized abelian surfaces over Fp and lifts to QQ

AUTHORS:

- Jean Kieffer (2025-06-17)

"""

# Copyright 2025 Jean Kieffer
# See LICENSE file for license details.

class PPASPrimeFiniteField(SageObject):
    r"""
    Create a data structure encoding a principally polarized abelian surface
    over a prime finite field Fp, and manipulatings lifts over QQ.
    """

    def __init__(self, data, inv_type = "modular"):
        r"""
        Initialize a PPASPrimeFiniteField object. The input data should be as
        specified by the PPASInvariants class, over a field GF(p) with p prime.
        """

        self.__fp_invariants = PPASInvariants.__init__(data, inv_type)
        self.__qq_lift = PPASRationalField.__init__([x.lift() for x in self.__fp_invariants.modular_invariants()])

    def isogenous_modular_invariants(self, ell, step = 1):
        pass
