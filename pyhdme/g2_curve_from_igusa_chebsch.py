
### Helper functions ###

# This is to prevent large error bounds when picking square roots of complex
# balls intersecting the real negative axis.
def g2_curve_safe_sqrt(x):
    try:
        if x.imag().contains_zero() and not x.real() >= 0:
            return I * sqrt(-x)
        else:
            return sqrt(x)
    except:
        return sqrt(x)

# This is to make sure quantities are actually nonzero even when manipulating
# complex balls
def g2_curve_safe_is_nonzero(x):
    try:
        return x.is_nonzero()
    except:
        return not x.is_zero()

### Detecting curves with specific automorphism groups ###

def g2_curve_has_generic_automorphisms(IC):
    R2 = g2_curve_R2_from_igusa_clebsch(IC)
    A, B, C, D = g2_curve_clebsch_from_igusa_clebsch(IC)
    return g2_curve_safe_is_nonzero(IC[3]) and g2_curve_safe_is_nonzero(R2) and (g2_curve_safe_is_nonzero(A)
                                                                                 or g2_curve_safe_is_nonzero(B)
                                                                                 or g2_curve_safe_is_nonzero(C))

def g2_curve_bolza_conditions_12(ABCD):
    return (ABCD[1] == 0) and (ABCD[2] == 0) and (ABCD[3] == 0)

def g2_curve_bolza_conditions_11(ABCD):
    t1 = 6 * ABCD[1] - ABCD[0]**2
    t2 = 6 * ABCD[2] - ABCD[0] * ABCD[1]
    return (t1 == 0) and (t2 == 0) and (ABCD[3] == 0)

def g2_curve_bolza_conditions_19(ABCD):
    t1 = - ABCD[1]**3 + 6 * ABCD[2]**2
    t4 = ABCD[0] * ABCD[1] + 6 * ABCD[2]
    t2 = -2 * t4 * ABCD[1] + 9 * ABCD[3]
    t3 = -15 * ABCD[2] + 2 * ABCD[0] * ABCD[1]
    r = (t1 == 0) and (t2 == 0)
    if r:
        if (ABCD[3] == 0 || t3 == 0):
            raise ValueError("Unexpected vanishing")
        a = 10 * g2_curve_safe_sqrt(-t4/t3)
    else:
        a = None
    return r, a

def g2_curve_bolza_conditions_23(ABCD):
    A, B, C, D = ABCD
    t1 = B**2 * A**3 - 6 * B * C + 4 * C * A**2 - 18 * D
    t2 = 4 * B**3 + 5 * C * B * A + 6 * C**2 - 3 * A * D
    t3 = 6 * C**2 - B**3
    r = (t1 == 0) and (t2 == 0)
    if r:
        if (D == 0 || t3 == 0):
            raise ValueError("Unexpected vanishing")
        a = (B**2 + A * C) / (2 * A**2 * B - 3 * A * C - 15 * B**2)
        a = 10 * g2_curve_safe_sqrt(a)
    else:
        a = None
    return r, a

### Cardona's reconstruction algorithm ###

def g2_curve_cardona_cubic(ABCD, conic, param):
    A, B, C, D = ABCD
    A11, A12, A22, A33 = conic

    a111 = 4/675*A^2*C - 8/225*B*C + 4/75*D
    a112 = 4/675*B^3 + 8/675*A*B*C + 8/225*C^2 + 2/225*A*D
    a122 = 2/675*A*B^3 + 8/2025*A^2*B*C + 8/675*B^2*C + 4/225*A*C^2 + 2/225*B*D
    a133 = -1/2025*A^2*B^4 - 4/6075*A^3*B^2*C + 8/2025*B^5 + 14/2025*A*B^3*C + 2/2025*A^2*B*C^2 + 8/675*B^2*C^2 + 4/675*A*C^3 + 1/225*A*B^2*D + 2/675*A^2*C*D + 2/225*B*C*D - 2/75*D^2
    a222 = 2/225*B^4 + 4/225*A*B^2*C + 16/2025*A^2*C^2 + 4/675*B*C^2 - 2/225*C*D
    a233 = 1/2025*A*B^5 + 2/1215*A^2*B^3*C + 8/6075*A^3*B*C^2 - 2/2025*B^4*C + 2/2025*A*B^2*C^2 + 8/2025*A^2*C^3 - 4/675*B*C^3 + 2/675*B^3*D + 1/675*A*B*C*D - 2/225*C^2*D - 1/225*A*D^2

    x, y, z = param
    t111 = -A33 * a111 * x**3
    t112 = -3 * A33 * a112 * x**2 * y
    t122 = -3 * A33 * a122 * x * y**2
    t133 = 3 * A22 * a133 * x * z**2
    t222 = -A33 * a222 * y**3
    t333 = 3 * A22 * a233 * x * z**2
    return t111 + t112 + t122 + t133 + t222 + t333

def g2_curve_from_igusa_clebsch_cardona(IC, t):
    A, B, C, D = g2_curve_clebsch_from_igusa_clebsch(IC)
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
    return g2_curve_cardona_cubic([A, B, C, D], [A11, A12, A22, A33], [x, y, z])

### Mestre's algorithm in the generic case ###

def g2_curve_mestre_conic(ABCD, U, I10):
    A, B, C, D = ABCD
    c11 = 2 * C + A * B/3
    c22 = D
    c33 = B * D/2 + 2 * C * (B**2 + A * C)/9
    c23 = B * (B**2 + A * C)/3 + C * (2 * C + A * B/3)/3
    c31 = D
    c12 = 2 * (B**2 + A * C)/3

    t11 = U**2 * I10**8 * c11
    t22 = I10**10 * c22
    t33 = U**8 * c33
    t23 = U**4 * I10**5 * c23
    t31 = U**5 * I10**4 * c31
    t12 = U * I10**9 * c12
    return [t11, t12, t13, t23, t31, t12]

def g2_curve_find_pt_on_conic(conic, t):
    for i in range(5**4):
        a2 = i % 5
        a3 = (i // 5) % 5
        b2 = 1 + (i // 25) % 5
        b3 = 1 + (i // 125) % 5

        x = t
        y = a2 * t + b2
        z = a3 * t + b3
        substitution = c11 * x**2 + c22 * y**2 + c33 * z**2 + c23 * y * z + c31 * x * z + c12 * x * y
        c0, c1, c2 = [substitution.coefficient(i) for i in range(3)]
        delta = c1**2 - 4 * c1 * c2
        if g2_curve_safe_is_nonzero(delta) and g2_curve_safe_is_nonzero(c2):
            delta = g2_curve_safe_sqrt(delta)
            x = (- c1 + delta) / (2 * c2)
            y = a2 * x + b2
            z = a3 * x + b3
            return [x, y, z]

    raise ValueError("Could not find point on conic")

def g2_curve_parametrize_conic(pt, conic, t):
    x0, y0, z0 = pt
    c11, c22, c23, c23, c31, c12 = conic

    #Enforce x0 != 0
    if not g2_curve_safe_is_nonzero(x0):
        if g2_curve_safe_is_nonzero(y0):
            y, z, x = g2_curve_parametrize_conic([y0, z0, x0], [c22, c33, c11, c31, c12, c23])
        elif g2_curve_safe_is_nonzero(z0):
            z, x, y = g2_curve_parametrize_conic([z0, x0, y0], [c33, c11, c22, c12, c23, c13])
        else:
            raise ValueError("Conic point does not have any nonzero coordinates")

    R = PolynomialRing(t.parent(), "u")
    x = x0
    y = y0 + u * t
    z = z0 + u
    substitution = c11 * x**2 + c22 * y**2 + c33 * z**2 + c23 * y * z + c31 * x * z + c12 * x * y
    a = substitution.coefficient(2) # in u
    b = substitution.coefficient(1) # in u
    return [x0 * a, y0 * a - t * b, z0 * a - b]

def g2_curve_mestre_cubic(ABCD, U, I10, param):
    A, B, C, D = ABCD
    c111 = 8 * (A**2 * C - 6 * B * C + 9 * D)/36
    c112 = 4 * (2 * B**3 + 4 * A * B * C + 12 * C**2 + 3 * A * D)/36
    c113 = 4 * (A * B**3 + 4 * A**2 * B * C/3 + 4 * B**2 * C + 6 * A * C**2 + 3 * B * D)/36
    c122 = 4 * (A * B**3 + 4 * A**2 * B * C/3 + 4 * B**2 * C + 6 * A * C**2 + 3 * B * D)/36
    c123 = 2 * (2 * B**4 + 4 * A * B**2 * C + 4 * A**2 * C**2/3 + 4 * B * C**2 + 3 * A * B * D + 12 * C * D)/36
    c133 = 2 * (A * B**4 + 4 * A**2 * B**2 * C/3 + 16 * B**3 * C/3 + 26 * A * B * C**2/3 +  8 * C**3 + 3 * B**2 * D + 2 * A * C * D)/36
    c222 = 4 * (3 * B**4 + 6 * A * B**2 * C + 8 * A**2 * C**2/3 + 2 * B * C**2 - 3 * C * D)/36
    c223 = 2 * (-2 * B**3 * C/3 - 4 * A * B * C**2/3 - 4 * C**3 + 9 * B**2 * D + 8 * A * C * D)/36
    c233 = 2 * (B**5 + 2 * A * B**3 * C + 8 * A**2 * B * C**2/9 + 2 * B**2 * C**2/3  - B * C * D + 9 * D**2)/36
    c333 = 1 * (-2 * B**4 * C - 4 * A * B**2 * C**2 - 16 * A**2 * C**3/9 - 4 * B * C**3/3  + 9 * B**3 * D + 12 * A * B * C * D + 20 * C**2 * D)/36

    x, y, z = param
    t111 = c111 * U**3 * I10**12 * x**3
    t112 = 3 * c112 * U**2 * I10**13 * x**2 * y
    t113 = 3 * c113 * U**6 * I10**8 * x**2 * z
    t122 = 3 * c122 * U * I10**14 * x * y**2
    t123 = 6 * c123 * U**5 * I10**9 * x * y * z
    t133 = 3 * c133 * U**9 * I10**4 * x * z**2
    t222 = c222 * I10**15 * y**3
    t223 = 3 * c223 * U**4 * I10**10 * y**2 * z
    t233 = 3 * c233 * U**8 * I10**5 * y * z**2
    t333 = c333 * U**12 * z**3
    return t111 + t112 + t113 + t122 + t123 + t133 + t222 + t223 + t233 + t333


def g2_curve_from_igusa_clebsch_mestre(IC, t):
    ABCD = g2_curve_clebsch_from_igusa_chebsch(IC)
    I2, I4, I6, I10 = IC
    A, B, C, D = ABCD
    if g2_curve_safe_is_nonzero(A):
        U = A**6
    elif g2_curve_safe_is_nonzero(B):
        U = B**3
    elif g2_curve_safe_is_nonzero(C):
        U = C**2
    else:
        raise ValueError("Could not find nonzero invariant of weight 12")

    conic = g2_curve_mestre_conic(ABCD, U, I10)
    pt = g2_curve_find_pt_on_conic(conic, t)
    x, y, z = g2_curve_parametrize_conic(pt, conic, t)
    return g2_curve_mestre_cubic(ABCD, U, I10, [x, y, z])

### main function ###

def g2_curve_from_igusa_clebsch(IC, base_ring=None, var_name="x"):

    if len(IC) != 4:
        raise ValueError("Input value should be 4 Igusa-Clebsch invariants")
    if IC[3] == 0:
        raise ValueError("Not the invariants of a genus 2 curve")

    R2 = g2_curve_R2_from_igusa_clebsch(IC)
    ABCD = g2_curve_clebsch_from_igusa_clebsch(IC)
    input_base_ring = IC[0].base_ring()
    try:
        input_is_exact = IC[0].is_exact() and IC[1].is_exact() and IC[2].is_exact() and IC[3].is_exact()
    except:
        input_is_exact = input_base_ring.is_exact()

    if base_ring is None:
        base_ring = input_base_ring
    poly_ring = PolynomialRing(base_ring, var_name)
    x = poly_ring.0

    if (g2_curve_has_generic_automorphisms(IC)):
        return g2_curve_from_igusa_clebsch_mestre(IC, x)
    elif not input_is_exact:
        raise ValueError("Cannot reconstruct curve with non-generic automorphisms on inexact input")
    elif g2_curve_safe_is_nonzero(R2):
        return x**6 + x
    elif g2_curve_bolza_conditions_12(ABCD):
        return x**5 + x
    elif g2_curve_bolza_conditions_11(ABCD):
        return x**6 + 1

    r, a = g2_curve_bolza_conditions_19(ABCD)
    if r:
        return x**6 + a * x**3 + 1
    r, a = g2_curve_bolza_conditions_23(ABCD)
    if r:
        return x**5 + a * x**3 + x

    return g2_curve_from_igusa_clebsch_cardona(IC, x)

