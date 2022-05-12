/* If flint was already previously included via another header (e.g.
 * flint_wrap.h) then it may be necessary to redefine ulong and slong again */

#ifndef ulong
#define ulong mp_limb_t
#define slong mp_limb_signed_t
#endif
