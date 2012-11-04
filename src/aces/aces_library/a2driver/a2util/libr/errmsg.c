
/*
 * prints the last error message
 */

#include <errno.h>
#include <stdio.h>
#include <string.h>

#include "f77_name.h"

void
F77_NAME(errmsg,ERRMSG)
()
{
    char *p;
/*
 * YAU - AIX 4.3 supports strerror just fine in string.h.
 * I don't know why this is here:
 *
 *  #ifdef _AIX
 *      perror("");
 *      p="strerror not supported on RS6000";
 *  #else
 *      p=strerror(errno);
 *  #endif
 *
 */
    p=strerror(errno);
    fprintf(stdout," System error: %d, %s\n",errno,p);
    fprintf(stderr," System error: %d, %s\n",errno,p);
    fflush(stdout);
    fflush(stderr);
}

