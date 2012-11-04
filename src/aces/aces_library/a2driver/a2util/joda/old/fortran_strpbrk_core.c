
/*
 * Fortran-callable routine to return the INDEX in sz of any single char
 * appearing in charset. Both sz and charset must be NULL-terminated.
 */

#include <string.h>

long fortran_strpbrk_core(const char * sz, const char * charset)
{
    char * cz = strpbrk(sz,charset);
    if (cz) return cz-sz+1;
    else    return 0;
}

long    fortran_strpbrk_core_(const char * sz, const char * charset)
{return fortran_strpbrk_core (             sz,              charset);}

