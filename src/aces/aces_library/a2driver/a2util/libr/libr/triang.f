      subroutine triang(ain,aout,nsiz,idiag)
C
C This routine accepts an NSIZ by NSIZ symmetric or antisymmetric matrix AIN
C and returns the its UPPER TRIANGLE in packed triangular form in AOUT.
C If IDIAG is set to 'Y', then the diagonal is included, if 'N' then 
C it is excluded. 
C
CEND
C      implicit double precision (a-h,o-z)
      implicit none
      integer nsiz
      double precision ain(nsiz,nsiz),aout(1)
      character *1 idiag
      character *(*) rtn
      parameter ( rtn = '  @uprtri')
C
      integer i, j, k
C
      integer ind, idum,jdum
      ind(idum,jdum)=(idum-k)*(idum-1)/2 + jdum
C
C      write (*,*) ' for i = 1, ', nsiz
      if (idiag .eq. 'N' .or. idiag .eq. 'n') then
C         write (*,*) ' No diagonals'
         k = 2
         do 40, i = 1, nsiz
            do 30, j = 1, i-1
C               write (*,*) i, j, ' | ', ind(i,j)
               aout(ind(i,j))=ain(j,i)
 30         continue
 40      continue
      elseif (idiag .eq. 'Y'  .or. idiag .eq. 'y') then
C         write (*,*)' Diagonals'
         k = 0
         do 20, i=1,nsiz
            do 10, j=1,i
C               write (*,*) i, j, ' | ', ind(i,j)
               aout(ind(i,j))=ain(j,i)
 10         continue
 20      continue
      else
         write (*,*) rtn,'-I: idiag must be ''Y'' or ''N'''
         write (*,*) rtn,'-F: Don''t recognize idiag ''',idiag,''''
         call errex
      endif
C
      return
      end
