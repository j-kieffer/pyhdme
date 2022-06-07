# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.


from sage.libs.flint.types cimport fmpz_poly_t, slong
from cysignals.signals cimport sig_on, sig_off
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz
from sage.libs.flint.fmpz_poly cimport fmpz_poly_init, fmpz_poly_clear, fmpz_poly_get_coeff_fmpz
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from cpython.int cimport PyInt_FromLong
from sage.libs.gmp.mpz cimport mpz_get_si
from sage.all import QQ, ZZ

def hecke_charpoly_wrapper(
    ell,
    wt,
    X
):
    r"""
    Given a prime number ell and an even weight wt >= 4, compute the characteristic polynomial of the Hecke operator of level ell on the space of Siegel modular forms of degree 2, full level, and weight wt, as a polynomial in X

    INPUT:

    - ``ell`` -- a prime number
    - ``wt`` -- an even number greater than 3
    - ``X`` -- a variable in a polynomial ring
    

    OUTPUT:

    - The characteristic polynomial of the Hecke operator, as a polynomial in X

    Examples::


    """
    
    ell = ZZ(ell)
    wt = ZZ(wt)
    assert ell.is_prime()
    
    cdef fmpz_poly_t charpoly
    cdef fmpz_t coeff
    cdef slong d
    cdef slong j
    
    cdef slong cell = mpz_get_si((<Integer>ell).value)
    cdef slong cwt = mpz_get_si((<Integer>wt).value)
    
    fmpz_poly_init(charpoly)
    fmpz_init(coeff)

    sig_on()
    assert hecke_charpoly(charpoly, cell, cwt) == 1
    sig_off()

    d = fmpz_poly_degree(charpoly)

    # convert charpoly to polynomial
    cdef object d_python
    d_python = PyInt_FromLong(d)
    
    ans = Rational(0) * X
    c = ZZ(0)

    for i in range(d_python+1):
        j = mpz_get_si((<Integer> ZZ(i)).value)
        fmpz_poly_get_coeff_fmpz(coeff, charpoly, j)
        fmpz_get_mpz((<Integer> c).value, coeff)
        ans += c * X**i

    fmpz_poly_clear(charpoly)
    fmpz_clear(coeff)
    return ans


