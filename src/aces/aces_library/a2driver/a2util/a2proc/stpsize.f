      Subroutine Stpsize(Omega, Deltaq, Delx, Nvibs)
C 
      Implicit Double Precision (A-H, O-Z)
C


c machsp.com : begin

c This data is used to measure byte-lengths and integer ratios of variables.

c iintln : the byte-length of a default integer
c ifltln : the byte-length of a double precision float
c iintfp : the number of integers in a double precision float
c ialone : the bitmask used to filter out the lowest fourth bits in an integer
c ibitwd : the number of bits in one-fourth of an integer

      integer         iintln, ifltln, iintfp, ialone, ibitwd
      common /machsp/ iintln, ifltln, iintfp, ialone, ibitwd
      save   /machsp/

c machsp.com : end









C
      Double Precision Na
C
      Dimension Deltaq(Nvibs), Omega(Nvibs)
     
      Parameter(h=6.62608D-34, c=2.99792458D+10, Na=6.02214D+23)
      
      PI = Acos(-1.0D+00)

      Do Imode = 1, Nvibs
         Deltaq(Imode) = Delx/(2.0D+00*PI*SQRT(c*Omega(Imode)/
     &                   h/Na*1.0D-23))
      Enddo 
C
      Return
      End

