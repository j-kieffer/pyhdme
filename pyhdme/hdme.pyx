# distutils: language=c
# clang c
# Copyright 2022 Edgar Costa
# See LICENSE file for license details.

from sage.libs.flint.types cimport fmpz, slong
from cysignals.signals cimport sig_on, sig_off
from memory_allocator cimport MemoryAllocator
from sage.libs.flint.fmpz cimport fmpz_init, fmpz_clear, fmpz_get_mpz, fmpz_set_mpz, fmpz_print
from sage.rings.integer cimport Integer
from sage.rings.rational cimport Rational
from cpython.int cimport PyInt_FromLong
from sage.libs.gmp.mpz cimport mpz_get_si
from sage.all import (
    ZZ,
    QQ,
    ceil,
    gcd,
    vector,
    prod,
)


from libc.stdio cimport printf
cdef print_vector(fmpz* elt, n):
    printf("(")
    for i in range(n):
        fmpz_print(&elt[i])
        if i != n-1:
            printf(", ");
    printf(")\n")

def rescale(c, I, weights):
    return vector([c**i * j for i, j in zip(weights, I)])

def canonicalize_rational_invariants(I, weights):
    assert len(I) == len(weights)
    I = [QQ(elt) for elt in I]
    weights = [ZZ(w) for w in weights] # for the divisions below
    # make it integral
    for i, w in enumerate(weights):
        if I[i] == 0:
            continue
        b, p = I[i].denominator().perfect_power()
        I = rescale(b**ceil(p/w), I, weights)
    assert all(elt.denominator() == 1 for elt in I)
    I = [ZZ(elt) for elt in I]
    d = gcd( j**(30/w) for w, j in zip(weights, I))
    c = prod( p**(-(e//30)) for p, e in d.factor())
    I = rescale(c, I, weights)
    oddw = [i for i, w in enumerate(weights) if w % 2 == 1]
    if len(oddw) > 0 and I[oddw[0]] < 0:
        I = rescale(-1, I, weights)
    return tuple(I)

def canonicalize_igusa_clebsch_invariants(I):
    return canonicalize_rational_invariants(I, [1,2,3,5])


def ic_from_igusa(M):
    M4, M6, M10, M12 = M
    I4 = M4*4
    I10 = -M10*2**12
    I12 = M12*2**15
    I2 = I12/I10
    I6prime = M6*4 # = (I2*I4-3*I6)/2
    I6 = (I2*I4 - I6prime*2)/3
    return canonicalize_igusa_clebsch_invariants((I2, I4, I6, I10))


def generic_wrapper(
    steps,
    igusa_clebsch_invariants,
    ell,
    verbose=False
):
    cdef int cverbose = int(verbose)
    for elt in [set_modeq_verbose, set_hecke_verbose]:
        elt(cverbose)
    assert steps in [1,2]
    assert len(igusa_clebsch_invariants) == 4
    ell = ZZ(ell)
    assert ell.is_prime()
    # https://beta.lmfdb.org/knowledge/show/g2c.igusa_clebsch_invariants
    igusa_clebsch_invariants = canonicalize_igusa_clebsch_invariants(igusa_clebsch_invariants)

    cdef MemoryAllocator mem = MemoryAllocator()
    cdef slong nb_roots
    cdef slong cell = mpz_get_si((<Integer>ell).value)

    # init output array
    max_nb_roots = (ell**3 + ell**2 + ell + 1)**steps
    cdef fmpz* all_isog_I = <fmpz*>mem.calloc(4 * max_nb_roots, sizeof(fmpz))
    for i in range(4 * max_nb_roots):
        fmpz_init(&all_isog_I[i])

    # convert igusa clebsch invariants to C
    cdef fmpz* cI = <fmpz*>mem.calloc(4, sizeof(fmpz))
    for i, elt in enumerate(igusa_clebsch_invariants):
        fmpz_init(&cI[i])
        fmpz_set_mpz(&cI[i], (<Integer>elt).value)
    # convert igusa clebsch invariants to igusa invariants
    # ie (I4/4, I6prime/4, -I10/2^12, I12/2^15)
    sig_on()
    igusa_from_IC_fmpz(cI, cI);
    sig_off()

    sig_on()
    if steps == 1:
        assert siegel_direct_isog_Q(&nb_roots, all_isog_I, cI, cell) == 1
    elif steps == 2:
        assert siegel_2step_direct_isog_Q(&nb_roots, all_isog_I, cI, cell) == 1
    sig_off()


    cdef object nb_roots_python
    nb_roots_python = PyInt_FromLong(nb_roots)
    assert nb_roots_python <= max_nb_roots;

    # convert all_isog_I to Igusa Clebsch
    # and convert it to python
    res = []
    sig_on()
    for i in range(nb_roots_python):
        #print_vector(&all_isog_I[4*i], 4)
        igusa_IC_fmpz(&all_isog_I[4*i], &all_isog_I[4*i])
        #print_vector(&all_isog_I[4*i], 4)
        #printf("##\n")
        r = [Integer(0) for _ in range(4)]
        for j in range(4):
            fmpz_get_mpz((<Integer>r[j]).value, &all_isog_I[i*4 + j])
        res.append(r)
    sig_off()

    for i in range(4):
        fmpz_clear(&cI[i])
    for i in range(4*max_nb_roots):
        fmpz_clear(&all_isog_I[i]);

    return [canonicalize_igusa_clebsch_invariants(elt) for elt in res]


def siegel_modeq_2step_isog_invariants_Q_wrapper(
      igusa_clebsch_invariants,
      ell,
      verbose=False
):
    r"""
    Given the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    of an abelian surface compute the Igusa-Clebsch invariants of surfaces (ell,ell,ell^2)-isogenous

    INPUT:

    - ``igusa_clebsch_invariants`` -- the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    - ``ell`` -- a prime number

    OUTPUT:

    - A list of Igusa-Clebsch invariants of abelian surfaces (ell, ell)-isogenous

    Examples::


    """
    return generic_wrapper(2, igusa_clebsch_invariants, ell, verbose=verbose)

def siegel_modeq_isog_invariants_Q_wrapper(
      igusa_clebsch_invariants,
      ell,
      verbose=False
):
    r"""
    Given the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    of an abelian surface compute the Igusa-Clebsch invariants of surfaces (ell,ell)-isogenous

    INPUT:

    - ``igusa_clebsch_invariants`` -- the Igusa-Clebsch invariants `I_2, I_4, I_6, I_{10}` of Igusa and Clebsch [IJ1960]_
    - ``ell`` -- a prime number

    OUTPUT:

    - A list of Igusa-Clebsch invariants of abelian surfaces (ell, ell)-isogenous

    Examples::


    """
    return generic_wrapper(1, igusa_clebsch_invariants, ell, verbose=verbose)

