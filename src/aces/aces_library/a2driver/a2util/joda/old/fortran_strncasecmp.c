
/*
 * Fortran-callable wrapper(s) for strncasecmp.
 */

#include <strings.h>

long fortran_strncasecmp (const char * s1, const char * s2, long * len)
{ return (long)strncasecmp(s1,s2,(size_t)*len); }

long fortran_strncasecmp_(const char * s1, const char * s2, long * len)
{ return (long)strncasecmp(s1,s2,(size_t)*len); }

