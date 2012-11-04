cvw   ===================================================================================
cvw
cvw   basic parameters:
cvw   maximum angular momentum of basis functions     (mxang)
cvw   maximum angular momentum of ECP projector + 1   (mxproj)
cvw
cvw   currently, up to i functions in basis   
cvw              up to h projectors          (note that this means 'lmax=5')
cvw
cvw   ===================================================================================
cvw      
      integer   mxang, mxproj
      parameter (mxang =6)
C if mxproj is to be changed, do remember to change the same variable in
C ecppar and ecpcab !!
      parameter (mxproj=5)
cvw      
cvw   ===================================================================================
cvw
cvw   derived parameters
cvw
cvw   lmaxcomb: maximum 'combined' angular momentum
cvw             either from basis function and ECP projector
cvw             or from the product of two basis functions.
cvw             This is the maximum angular momentum for which
cvw             angular momentum gymnastics have to be performed
      integer lmaxcomb
      parameter(lmaxcomb=mxang+max(mxang,mxproj-1))
cvw
cvw   nftmax:   maximum number of cartesian mononomials
cvw             (i. e. the dimension of the n_cart, l_cart, m_cart arrays)
cvw   nftmxg:   the same for gradient calculations
cvw
      integer nftmax,nftmxg
cvw
cvw   ndegen(l) = (l+1)(l+2)/2, number of cartesian components for ang. mom. l
cvw   nftmax = Sum{l=0 ... mxang} ndegen(l)
cvw   nftmxg = Sum{l=0 ... mxang} ndegen(l-1)+ndegen(l+1)
cvw   nftmxe = Sum{l=0 ... mxang} ndegen(l+2)
cvw
cvw   nftmxg primitives arise in geo gradient
cvw   nftmxe primitives arise in basis gradient (exponent derivatives)
cvw
      parameter (nftmax=(mxang+1)*(mxang+2)*(mxang+3)/6)
      parameter (nftmxg=(mxang+1)*(9+5*mxang+mxang*mxang)/3)
      parameter (nftmxe=(mxang+1)*(36+11*mxang+mxang*mxang)/6)
cvw 
cvw   ===================================================================================
cvw 
cvw   TABLES FOR ANGULAR MOMENTUM GYMNASTICS:
cvw
cvw   angular momentum data is indexed by l, m
cvw   l=1,2,3,.... ,lfdim
cvw   m = 1,2,3,... 2*l-1
cvw
cvw   and stored in a linear array at position lf(l)+m
cvw   linear dimension is lmfdim = lfdim**2
cvw
cvw   Real spherical harmonics Y(l,m) are expressed as a linear combination
cvw   of cartesian products ZLM(i,j,k) x**i y**j z**k
cvw
cvw   The coefficients are in zlm(istart ... iend)
cvw   with istart=lmf(lf(l)+m), iend=lml(lf(l)+,m)
cvw   and the corresponding values of i,j,k in the arrays lmx, lmy, lmz
cvw
cvw   Note that for a given l,m, all cartesian terms have the same
cvw   parity pattern of the ijk. This pattern is stored in lmm(lf(l)+m)
cvw
cvw   bit 0 of lmm is set if k is EVEN
cvw   bit 1 of lmm is set if j is EVEN
cvw   bit 2 of lmm is set if i is EVEN
cvw
cvw   The total length of the zlm, lmx, lmy, lmz arrays is lmxdim, a number
cvw   which must be calculated by a rather complex formula
cvw
cvw   There are some auxiliary tables such as double factorial (dfac),
cvw   its inverse (dfaci) and binomial coefficients (binom_coef)
cvw   which are also stored in this data structure.
cvw
      integer lfdim,lmfdim,lmxdim
      parameter(lfdim=lmaxcomb+1)
      parameter(lmfdim=lfdim**2)
      parameter(lmxdim=(
     &    lmaxcomb*(lmaxcomb+2)*(lmaxcomb+3)*(lmaxcomb+4)/3 +
     &    (lmaxcomb+2)**2 * (lmaxcomb+4)
     &    ) /16)
cvw      
cvw   ===================================================================================
cvw
cvw   cvw_ecp_r/cvw_ecp_i: readonly data, loaded by subr. tab_ecp
cvw                        note that flmtx, mc, mr, mrclo, mrchi
cvw                        are only needed for spin-orbit integrals
cvw                        (compressed matrices of Lx, Ly, Lz between
cvw                         real spherical harmonics of same l)
cvw
cvw   common cvw_ecp_i contains the integer/logical data
cvw   common cvw_ecp_r contains the floating point data
cvw
      double precision zlm(lmxdim)
      double precision flmtx(mxproj*mxproj,3)
      double precision binom_coef(lfdim*(lfdim+1)/2)
      integer          lf(lfdim),lmf(lmfdim),lml(lmfdim)
      integer          lmx(lmxdim),lmy(lmxdim),lmz(lmxdim)
      integer          lmm(lmfdim)
      integer          n_cart(nftmax),l_cart(nftmax),m_cart(nftmax)
      integer          mc(mxproj*mxproj,3),mr(mxproj*mxproj,3)
      integer          mrclo(mxproj),mrchi(mxproj)
      double precision dfac (4*mxang+2*mxproj+9)
      double precision dfaci(4*mxang+2*mxproj+9)
      double precision fprod(lfdim,lfdim)
cvw      
      integer          nprim_shell(mxang+1)
      integer          nprim_shellg(mxang+1)
      integer          nlm_grd(nftmxg,7)
      integer          nlm_min(mxang+1)
      integer          nlm_max(mxang+1)
      integer          nlm_ming(mxang+1)
      integer          nlm_maxg(mxang+1)
cvw
cvw   mrclo, mrchi, mc, mr:    only used for SO integrals
cvw   nlm_grd:                 only used in gradient code
cvw
      common /cvw_ecp_r/ zlm,binom_coef,flmtx,fprod,
     &                   dfac,dfaci
      common /cvw_ecp_i/ lf,lmf,lml,lmx,lmy,lmz,lmm,
     &                   n_cart,l_cart,m_cart,
     &                   mrclo, mrchi, mc, mr,
     &                   nlm_grd,nlm_min,nlm_max,
     &                   nlm_ming,nlm_maxg,
     &                   nprim_shell,nprim_shellg
cvw      
cvw   ===================================================================================
