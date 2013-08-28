      subroutine triangs(ain,aout,nsiz,idiag,half)
C
C This routine accepts an NSIZ by NSIZ symmetric or antisymmetric matrix AIN
C and returns its upper or lower triangle  in packed triangular form in AOUT.
C
C IDIAG   'Y'  diagonal is included.
C         'N'  diagonal is excluded.
C
C HALF    'U'  upper triangle is returned.
C         'L'  lower triangle is returned.
C
CEND
C RPM March 1994
C      implicit double precision (a-h,o-z)
      implicit none
      integer nsiz
      double precision ain(nsiz,nsiz),aout(1)
      character *1 idiag,half
      character *(*) rtn
      parameter ( rtn = '  @lwrtri')
C
      integer i, j, k
C
      integer ind, idum,jdum
      ind(idum,jdum)=(idum-k)*(idum-1)/2 + jdum
C
C      write (*,*) ' for i = 1, ', nsiz
      if ( half .eq. 'L'  .or.  half .eq. 'l' ) then
C         write (*,*) ' lower triangle'
         if (idiag .eq. 'N' .or. idiag .eq. 'n') then
C            write (*,*) ' No diagonals'
            k = 2
            do 20, i = 1, nsiz
               do 10, j = 1, i-1
C               write (*,*) i, j, ' | ', ind(i,j)
                  aout(ind(i,j))=ain(i,j)
 10            continue
 20         continue
         elseif (idiag .eq. 'Y' .or. idiag .eq. 'y') then
C         write (*,*)' Diagonals'
            k = 0
            do 40, i=1,nsiz
               do 30, j=1,i
C               write (*,*) i, j, ' | ', ind(i,j)
                  aout(ind(i,j))=ain(i,j)
 30            continue
 40         continue
         else
            write (*,*) rtn,'-I: idiag must be ''Y'' or ''N'''
            write (*,*) rtn,'-F: Don''t recognize idiag ''',idiag,''''
            call errex
         endif
      elseif ( half .eq. 'U' .or. half .eq. 'u' ) then
C         write(*,*) ' Upper triangle'
         if (idiag .eq. 'N' .or. idiag .eq. 'n') then
C         write (*,*) ' No diagonals'
            k = 2
            do 60, i = 1, nsiz
               do 50, j = 1, i-1
C               write (*,*) i, j, ' | ', ind(i,j)
                  aout(ind(i,j))=ain(j,i)
 50            continue
 60         continue
         elseif (idiag .eq. 'Y'  .or. idiag .eq. 'y') then
C         write (*,*)' Diagonals'
            k = 0
            do 80, i=1,nsiz
               do 70, j=1,i
C               write (*,*) i, j, ' | ', ind(i,j)
                  aout(ind(i,j))=ain(j,i)
 70            continue
 80         continue
         else
            write (*,*) rtn,'-I: IDIAG must be ''Y'' or ''N'''
            write (*,*) rtn,'-F: Don''t recognize IDIAG ''',idiag,''''
            call errex
         endif
      else
         write (*,*) rtn,'-I: HALF must be ''U'' or ''L'''
         write (*,*) rtn,'-F: Don''t recognize HALF ''',half,''''
         call errex
      endif
C
      return
      end
