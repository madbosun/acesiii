
#ifdef __fortran77

C THIS ROUTINE COPIES THE FIRST LEN ELEMENTS OF INTEGER
C VECTOR I1 INTO VECTOR I2. IT IS DESTRUCTIVE IF I1 AND I2 OVERLAP.

      SUBROUTINE ICOPY(I1,I2,LEN)
      INTEGER I1(LEN),I2(LEN)
      DO I=1,LEN
         I2(I)=I1(I)
      END DO
      RETURN
      END

#endif /* __fortran77 */

#ifdef __ansi_c
#include <string.h>
void
#ifdef _UNICOS
     ICOPY
#else
#  ifdef C_SUFFIX
     icopy_
#  else
     icopy
#  endif /* C_SUFFIX */
#endif /* _UNICOS */
(unsigned char * s1, unsigned char * s2, long * length)
{
    long i;
    memcpy(s2,s1,(size_t)(*length*sizeof(i)));
    return;
}
#endif /* __ansi_c */

