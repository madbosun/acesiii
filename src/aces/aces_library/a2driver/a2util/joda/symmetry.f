













































































































































































































































































































































































































































































































































      SUBROUTINE SYMMETRY(SCRATCH,QTMP,NEWQ,NOSILENT)
C
C This is a front end routine that handles geometry/symmetry related
C issues. With this, we can smoothly turn off the symmetry
C processing when it is requested. Also, the logic that pertains
C to the reorientation is much simpler. Ajith Perera, 11/05
C 
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C MXATMS     : Maximum number of atoms currently allowed
C MAXCNTVS   : Maximum number of connectivites per center
C MAXREDUNCO : Maximum number of redundant coordinates.
C
      INTEGER MXATMS, MAXCNTVS, MAXREDUNCO
      PARAMETER (MXATMS=200, MAXCNTVS = 10, MAXREDUNCO = 3*MXATMS)
c     Maximum string length of absolute file names
      INTEGER FNAMELEN
      PARAMETER (FNAMELEN=80)
C coord.com : begin
C
      DOUBLE PRECISION Q, R, ATMASS
      INTEGER NCON, NR, ISQUASH, IATNUM, IUNIQUE, NEQ, IEQUIV,
     &        NOPTI, NATOMS
      COMMON /COORD/ Q(3*MXATMS), R(MAXREDUNCO), NCON(MAXREDUNCO),
     &     NR(MXATMS),ISQUASH(MAXREDUNCO),IATNUM(MXATMS),
     &     ATMASS(MXATMS),IUNIQUE(MAXREDUNCO),NEQ(MAXREDUNCO),
     &     IEQUIV(MAXREDUNCO,MAXREDUNCO),
     &     NOPTI(MAXREDUNCO), NATOMS

C coord.com : end


C cbchar.com : begin
C
      CHARACTER*5 ZSYM, VARNAM, PARNAM
      COMMON /CBCHAR/ ZSYM(MXATMS), VARNAM(MAXREDUNCO),
     &                PARNAM(MAXREDUNCO)

C cbchar.com : end



      DOUBLE PRECISION IT(3,3),CM(3),IV(3,3),NEWQ(NX),dtmp
      DOUBLE PRECISION RM(3,3),QTMP(NX),TATB(3),QOLD(100)
      DOUBLE PRECISION SCRATCH(3*NX),MOLWT,ORIEN2(9)
      LOGICAL XYZIN,NWFINDIF,NOSILENT
      INTEGER IORGRP, MEMBER(MXATMS)
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /TOLERS/ SYMTOL,DEGTOL
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $     ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $     XYZTol
      COMMON /INPTYP/ XYZIN, NWFINDIF
      Character*(8*mxatms) szStSymTmp
      Character*8 STSYM(MXATMS)
C
      CHARACTER*(5*MXATMS) ZSYMUNI
C
c io_units.par : begin

      integer    LuOut
      parameter (LuOut = 6)

      integer    LuErr
      parameter (LuErr = 6)

      integer    LuBasL
      parameter (LuBasL = 1)
      character*(*) BasFil
      parameter    (BasFil = 'BASINF')

      integer    LuVMol
      parameter (LuVMol = 3)
      character*(*) MolFil
      parameter    (MolFil = 'MOL')
      integer    LuAbi
      parameter (LuAbi = 3)
      character*(*) AbiFil
      parameter    (AbiFil = 'INP')
      integer    LuCad
      parameter (LuCad = 3)
      character*(*) CadFil
      parameter    (CadFil = 'CAD')

      integer    LuZ
      parameter (LuZ = 4)
      character*(*) ZFil
      parameter    (ZFil = 'ZMAT')

      integer    LuGrd
      parameter (LuGrd = 7)
      character*(*) GrdFil
      parameter    (GrdFil = 'GRD')

      integer    LuHsn
      parameter (LuHsn = 8)
      character*(*) HsnFil
      parameter    (HsnFil = 'FCM')

      integer    LuFrq
      parameter (LuFrq = 78)
      character*(*) FrqFil
      parameter    (FrqFil = 'FRQARC')

      integer    LuDone
      parameter (LuDone = 80)
      character*(*) DonFil
      parameter    (DonFil = 'JODADONE')

      integer    LuNucD
      parameter (LuNucD = 81)
      character*(*) NDFil
      parameter    (NDFil = 'NUCDIP')

      integer LuFiles
      parameter (LuFiles = 90)

c io_units.par : end
C     Symmetry Information
C     FPGrp   Full point group
C     BPGrp   Largest Abelian subgroup
C     PGrp    "Computational" point group
      Character*4 FPGrp, BPGrp, PGrp, PTGrp
      Character*8 TMPGrp
      Common /PtGp_com/ FPGrp, BPGrp, PGrp
      Common /Orient/ Orient(3,3)


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



      DATA IONE, ININE /1, 9/
      DATA ONE   / 1.0D+00/
      DATA ONEM  /-1.0D+00/
      DATA ZILCH / 0.0D+00/
      DTOR=DACOS(ONEM)/180.D0
C
C TRANSLATE TO CENTER OF MASS
C
      IDEGEN=0
C
      IF(IPRNT.GE.3 .AND. NOSILENT) WRITE(LUOUT,7733)
     &                              (ATMASS(J),J = 1,NATOMS)
 7733 FORMAT(1X,F15.10)
C
      Write(6,"(a, 5I)"), "The number of atoms", Natoms
      Write(6,*)
      Print*, "The Cartesians before translations to CM" 
      Write(6, "(3F10.5)"), (Q(I), I=1, 3*NATOMS)
      Write(6,*)
      CALL DCOPY(3*NATOMS, Q, 1, QOLD, 1)
C
      CMX=0.D0
      CMY=0.D0
      CMZ=0.D0
      MOLWT=0.D0
      DO 20 I = 1,NATOMS
        CMX = ATMASS(I)*Q(3*I-2)+CMX
        CMY = ATMASS(I)*Q(3*I-1)+CMY
        CMZ = ATMASS(I)*Q(3*I)+CMZ
   20 MOLWT = MOLWT+ATMASS(I)
      IF (MOLWT.LT.1.0D-8) THEN
         WRITE(LUOUT,*) '@SYMMETRY: No real atoms in Z-matrix.'
cYAU         CALL GFNAME(ARCFIL,FNAME,ILENGTH)
cYAU         OPEN(UNIT=LUARC,File=fname(1:ilength),STATUS='OLD')
cYAU         CLOSE(UNIT=LUARC,STATUS='DELETE')
         CALL ERREX
      END IF
      CM(1) = CMX/MOLWT
      CM(2) = CMY/MOLWT
      CM(3) = CMZ/MOLWT
      DO 30 I = 1,NATOMS
        DO 301 J = 0,2
          Q(3*I-J) = Q(3*I-J)-CM(3-J)
  301   continue
   30 continue
      Write(6,*), "The Cartesians in center of mass coords"
      Write(6, "(3F12.7)"), (Q(I), I=1, 3*NATOMS)
      Write(6,*)
      IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,*)
     &     'After translation to center of mass coordinates '
      IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,80)(Q(I),I = 1,NX)
   80 FORMAT((4X,3(2X,F16.12)))
c
cjdw 5/26/95
c
      IF (IPRNT .GE. 4 .AND. NOSILENT) THEN
      write(6,*)
      write(6,*) ' @symmetry-i, Coordinates after  COM shift '
      do 1010 i=1,natoms
        write(6,'(3F20.12)') q(3*i-2),q(3*i-1),q(3*i)
 1010 continue
      write(6,*)
      ENDIF
C
C     BUILD INERTIA TENSOR
C
      CALL INERTIA(IT, NOSILENT)
      IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,*)'Inertia tensor:'
      IF(IPRNT .GE. 4 .AND. NOSILENT)WRITE(LUOUT,80)((IT(I,J),J = 1,3),
     &                               I = 1,3)
C
C DIAGONALIZE INERTIA TENSOR.
C
      CALL FILTER(IT,9,1.D-12)
      CALL EIG(IT,IV,3,3,0)
      CALL SCOPY(9,IV,1,ORIEN2,1)
C
C CHECK *NOW* FOR DEGENERACY OF EIGENVALUES -- IF PRESENT, THEN SEE IF
C  UNIQUE MOMENT OF INERTIA IS ALONG X.  IF SO, CHANGE THIS AXIS TO Z
C  BY ROTATING THE EIGENVECTOR MATRIX ABOUT Y.  THIS GUARANTEES THAT
C  THE UNIQUE AXIS WILL LIE ALONG Z.
C
      ABSDIF=ABS(IT(2,2)-IT(3,3))
      Z=ABS(MAX(IT(2,2),IT(3,3)))
      IF(Z.NE.0.0)THEN
        RELDIF=ABSDIF/Z
      ELSE
        RELDIF=0.0
      ENDIF
      IF(RELDIF.LT.DEGTOL)THEN
        RANG=90.D0
        CALL ROTM(2,RANG,1,RM)
        CALL MODMATMUL(IV,IV,RM,3,3,3,3,3,3)
        ATMP=IT(1,1)
        IT(1,1)=IT(3,3)
        IT(3,3)=ATMP
      ENDIF

CJDW 1/6/98. Replace MATMULV by XGEMM call. NEWQ = (IV)^T * Q.
C     CALL MATMULV(NEWQ,Q,IV,NATOMS,3,3)
      CALL XGEMM('T','N',3,NATOMS,3,ONE,IV,3,Q,3,ZILCH,NEWQ,3)
      IF(IPRNT .GE. 4 .AND. NOSILENT)THEN
        WRITE(LUOUT,*)'   Diagonalized inertia tensor:'
        WRITE(LUOUT,80)((IT(I,J),J = 1,3),I = 1,3)
        WRITE(LUOUT,*)'   Eigenvectors of inertia tensor: '
        WRITE(LUOUT,80)((Iv(I,J),J = 1,3),I = 1,3)
        WRITE (LUOUT,*) '   Principal axis orientation ',
     $       'for molecular system: '
        WRITE (LUOUT,80) (NEWQ(I),I = 1,NX)
      ENDIF
C
C     Print out the rotational constants from the inertia
C
      Call RotCon (IT, 1, IErr, NOSILENT)
C
C CHECK HANDEDNESS OF INERTIAL AXES AND SWAP IF NECESSARY.  DETERMINANT
C   HAS TO BE UNITY.
C
      CALL CROSS(IV(1,1),IV(1,2),TATB,1)
      XXX=xdot(3,TATB,1,IV(1,3),1)
      IF(IPRNT.GE.4 .AND. NOSILENT)THEN
        WRITE(LUOUT,3143)XXX
 3143   FORMAT(T3,'@SYMMETRY-I, Handedness of inertial frame:',F8.5)
      ENDIF
      IF(XXX.LT.0.D0)THEN
        IF(IPRNT.GE.4 .AND. NOSILENT)WRITE(LUOUT,3142)
 3142   FORMAT(T3,'@SYMMETRY-I, Sense of inertial frame ',
     &       'y-axis will be changed.')
        DO 1451 I=2,NX-1,3
 1451   NEWQ(I)=-NEWQ(I)
      ENDIF
C
C First, check the degeneracies if the inertia matrix to identify 
C whether we have a degenracy in the inertia matrix. When there
C is no degenaracy, we know immediately that the molecule has 
C abelian point group (see symmetry_auto for symmetry processing). 
C Also, we can assess the linearity by looking at the diagonal of
C the moments of inertia matrix.
C
         ILINEAR = 0
         NP1  = NATOMS+1
         JAX  = NX+1
         JAX2 = NX+JAX
         X    = ONEM
         DO I = 1, 3
            Z      = IT(I,I)
            ABSDIF = ABS(Z-X)
            Z0     = ABS(MAX(X,Z))
            IF (Z0.LT.1.0D-14) THEN
               IF (ABSDIF.LT.1.0D-14) THEN
                  RELDIF=0.0D0
               ELSE
                  RELDIF=1.0D0
               END IF
            ELSE
               RELDIF=ABSDIF/Z0
            END IF
            IF (RELDIF.LT.DEGTOL) IDEGEN = IDEGEN + 1
            X = Z
            IF (DABS(IT(I,I)).LT.SYMTOL) ILINEAR=1
         END DO

         IF (IPRNT.GE.4 .AND. NOSILENT) WRITE(6,3146) IDEGEN+1
 3146    FORMAT(T3,'@SYMMETRY-I, The symmetry group is ',i1,'-fold ',
     &          'degenerate.')
C
      Write(6,*)
      Print*, "The symmetry processing begins"

      IF (IFLAGS(60).EQ.0) THEN
      Print*, "The symmetry is not used for anything"
      Write(6,*)
C        
C Undo the rotation and translation
C 
CSSS      CALL XGEMM('N','N',3,NATOMS,3,ONE,IV,3,NEWQ,3,ZILCH,Q,3)
CSSS      DO I = 1,NATOMS
CSSS         DO J = 0,2
CSSS            Q(3*I-J) = Q(3*I-J)+ CM(3-J)
CSSS         ENDDO
CSSS      ENDDO
      Print*, "The Cartesians wrt CM origin"
      Write(6, "(3F10.5)"), (Q(I), I=1, 3*NATOMS)
      Write(6,*)
c
C The user does not want to use symmetry. Therefore, skip
C all symmetry-related work beside automatically setting the
C full point group, computational point, and largest abelian subgroup.
C (all are C1). Let's re-orient the molecule to the principal
C axis coordinate system unless the user specificaly says NO.
C
         FPGRP = "C1  "
         PGRP  = "C1  "
         BPGRP = "C1  "
         DO I = 1, 3
            ORIENT(I,I) = ONE
         END DO
         CALL DCOPY(9, ORIEN2, 1, SCRATCH, 1)
C
C Since under this condition no call are made to SYMDRV.F, so all
C the JOBACR records that are created by FLUSH.F (call by SYMDRV.F)
C recorded. The content of these records can be trivialy generated
C for SYM=NONE (same as SYM=C1). 
C 
         CALL PUTREC(20,'JOBARC','LINEAR  ',IONE,ILINEAR)
         CALL PUTREC(20,'JOBARC','FULLNORB',IONE,NATOMS)
         CALL PUTREC(20,'JOBARC','COMPNORB',IONE,NATOMS)
C
         ISIZE=NATOMS
         iNdx = IONE
         do i = 1, iSize
             STSYM(i)  = 'C1      '
             MEMBER(i) = i
             szStSymTmp(iNdx:iNdx+7) = STSYM(i)(1:8)
           iNdx = iNdx + 8
         end do
         CALL PUTCREC(20,'JOBARC','FULLSTGP',iSize*8,
     &                szStSymTmp(1:iSize))
         CALL PUTCREC(20, 'JOBARC', 'FULLPTGP', 8, STSYM(1))
         CALL PUTCREC(20, 'JOBARC', 'COMPPTGP', 8, STSYM(1))
         CALL PUTREC(20,'JOBARC','COMPMEMB',ISIZE,MEMBER)
         CALL PUTREC(20,'JOBARC','FULLMEMB',ISIZE,MEMBER)
         CALL PUTREC(20,'JOBARC','COMPPERM',ISIZE,MEMBER)
         CALL PUTREC(20,'JOBARC','FULLPERM',ISIZE,MEMBER)
         DO i = 1, iSize
            MEMBER(i) = 1
         ENDDO
         CALL PUTREC(20,'JOBARC','COMPPOPV',ISIZE,MEMBER) 
         CALL PUTREC(20,'JOBARC','FULLPOPV',ISIZE,MEMBER) 
         CALL PUTREC(20,'JOBARC','COMPNIRR',IONE,IONE)
         CALL PUTREC(20,'JOBARC','FULLNIRR',IONE,IONE)
         CALL PUTREC(20,'JOBARC','FULLORDR',IONE,IONE)
         CALL PUTREC(20,'JOBARC','COMPORDR',IONE,IONE)
         CALL IDNMAT(RM,3,IJUNK)
         CALL PUTREC(20,'JOBARC','FULLSYOP',ININE*IINTFP,RM) 
         CALL PUTREC(20,'JOBARC','COMPSYOP',ININE*IINTFP,RM) 
C
         CALL GETREC(-1, 'JOBARC', '12SWITCH', ione, ibad)
         II = 0
         DO I = 1, NATOMS
C
C Minor fix to get Symmetry=none to work with dummy atoms; 02/2012
C Ajith Perera
C
            IF (IATNUM(I) .EQ. 0) Then
               MEMBER(I) = 999
            ELSE
               II = II + 1
               MEMBER(I) = II
            ENDIF
         ENDDO

         IF (IBAD .EQ. 1) THEN 
            MEMBER(2) = 1
            MEMBER(1) = 2 
         ENDIF
         CALL PUTREC(20,'JOBARC','ZMAT2MOL',ISIZE,MEMBER)
C
C Copy the principal axis orientation (NEWQ to Q) and print it.
C It will be put into JOBARC (see below). 
C
         IF (IFLAGS2(125).EQ.1) THEN
C
C Undo the rotation and translation, if xxx less than zero, then 
C Handedness of the inertial frame has changed and we must undo
C that before undo the principal axix rotation (See above). 
C
            IF (XXX .LT. 0.0D0) Then
               DO I=2,NX-1,3
                  NEWQ(I)=-NEWQ(I) 
               ENDDO
            ENDIF  
C
            CALL XGEMM('N','N',3,NATOMS,3,ONE,IV,3,NEWQ,3,ZILCH,Q,3)
      Write(6,*)
      Print*, "The Cartesians after no-principal axis rot. "
      Write(6, "(3F10.5)"), (Q(I), I=1, 3*NATOMS)
            DO I = 1,NATOMS
               DO J = 0,2
                  Q(3*I-J) = Q(3*I-J)+ CM(3-J)
               ENDDO
            ENDDO
      Write(6,*)
      Print*, "The Cartesians after no-reorintation CM"
      Write(6, "(3F10.5)"), (Q(I), I=1, 3*NATOMS)
         ELSE
C
C When symmetry is none NOREORI is turned on by default and so
C in principle this block is inactive.
C 
            CALL DCOPY(3*nAtoms, NEWQ, 1, Q, 1)

      Write(6,*)
      Print*, "The Cartesians used when sym=none,noreori=off"
      Write(6, "(3F10.5)"), (NEWQ(I), I=1, 3*NATOMS)
         END IF
C
         IF (IPRNT .GT. 4 .AND. NOSILENT) THEN
             WRITE(6,1040) 
 1040        FORMAT(/,' @SYMMETRY-I, ',
     &              'Cartesian coordinates (Bohr): ',/)
             do i = 1, natoms
                write(6,'(3F20.12)') q(3*i-2),q(3*i-1),q(3*i)
             end do
         ENDIF
        IF (NOSILENT) THEN
           WRITE(LUOUT,788)
           WRITE(LUOUT,789)
  789      FORMAT('The symmetry processing is fully turned',
     &            ' off and no symmetry opertion are',
     &            ' applied.')
           WRITE(LUOUT, *)
           WRITE(LUOUT,1881)FPGRP
 1881      FORMAT(T3,' The full molecular point group is ',a,'.')
           WRITE(LUOUT,177)BPGRP
           WRITE(LUOUT,712)PGRP
           WRITE(LUOUT,788)
           WRITE(LUOUT,*)
  788      FORMAT(80('*'))
  712      FORMAT(T3,' The computational point group is ',A,'.')
  177      FORMAT(T3,' The largest Abelian subgroup of the full ',
     $               'molecular point group is ',a,'.')
        ENDIF
C
      ELSE
C
      Print*, "The Cartesians before entering symmetry auto"
      Write(6, "(3F10.5)"), (NEWQ(I), I=1, 3*NATOMS)
      Write(6, "(3F10.5)"), (QTMP(I), I=1, 3*NATOMS)
 

            CALL SYMMETRY_AUTO(SCRATCH, QTMP, NEWQ, IT, IDEGEN, 
     &                         ORIEN2, NOSILENT)
C
C All of the no reorientation stuff should be here. Two things can happen:
C molecule can be in principal axis orientation, or it could be in
C another orientation dictated by the point group of the molecule.
C
         IF (IFLAGS2(125).EQ.1) THEN
            IF (PGRP.EQ."C1  ") THEN
C
C Undo the rotation and translation, if xxx less than zero, then
C Handedness of the inertial frame has changed and we must undo
C that before undo the principal axix rotation (See above).
C
            IF (XXX .LT. 0.0D0) Then
               DO I=2,NX-1,3
                  NEWQ(I)=-NEWQ(I)
               ENDDO
            ENDIF

      Write(6,*), "The Cartesians in the principal axis rotation ori."
      Write(6, "(3F12.7)"), (NEWQ(I), I=1, 3*NATOMS)                  
      Write(6,*)
      Write(6,*) "The transformation matrix"
      Write(6,"(3F12.7)")((Iv(I,J),J = 1,3),I = 1,3)
      Write(6,*)
C
C Reorientation is due to rotation of the principal axis.
C If no reorientation is requested, undo the rotation and translation. 
C
               CALL XGEMM('N','N',3,NATOMS,3,ONE,IV,3,NEWQ,3,ZILCH,
     &                     QTMP,3)
      Write(6,*), "The Cartesians in center of mass coords"
      Write(6, "(3F12.7)"), (QTMP(J), J=1, 3*NATOMS)                  
      Write(6,*)
               DO I = 1,NATOMS
                  DO J = 0,2
                     QTMP(3*I-J) = QTMP(3*I-J)+CM(3-J)
                  ENDDO
               ENDDO
               CALL PUTREC(20,'JOBARC','NOREOCOR',3*NATOMS*IINTFP,
     &                     QTMP)
               CALL DCOPY(3*NATOMS, QTMP, 1, Q, 1)

      Write(6,*), "The Cartesians after undoing the sym. REORIENT"
      Write(6, "(3F12.7)"), (QTMP(J), J=1, 3*NATOMS)                  
C
            ELSE
C
C We can stop the reorientation (perhaps easily), but let's leave that
C for another round. I need to know something about the irrep
C labels in the rest of the code. Best thing to do is to do the
C calculation in the original geometry but leave the symmetry
C labeling to match with the rotated geometry.
C
               CALL PUTREC(20,'JOBARC','NOREOCOR',3*NATOMS*IINTFP,
     &                     QOLD)

      Print*, "The Cartesians before translations to CM:NO REORIENT"
      Print*, (QOLD(I), I=1, 3*NATOMS)
C
            END IF
         END IF
C
C No reorientation is requested and the symmetry can be C1 or any
C othe point group. For those cases let it be that it always use
C center of mass/principal axis defined with respect to the moment of inertia 
C tensor as the base coordinate system. When there is symmetry, the
C molecule is oriented according to the standard group theory conventions.
C THE WORKING COORDINATES ARE IN THE ARRAY Q of the /COORD/ common block.
C NOTE THAT THE COORDS IN Q ARE ORIENTED ACCORDING TO GROUP THOERY CONVENTIONS
C AND THE ACTUAL PROCESS OF REORIENTING HAPPENS IN SYMMETRY_AUTO AND 
C SYMMETRY_AUTO IS CALLED FOR BOTH SYMMETRY{ON, OFF}.
C
      END IF
C
C DUMP SOME NECESSARIES TO JOBARC
C
      CALL PUTREC(20,'JOBARC','ORIENT2 ',ININE*IINTFP,SCRATCH)
      CALL PUTREC(20,'JOBARC','NATOMS  ',IONE,NATOMS)
      CALL PUTREC(20,'JOBARC','COORD   ',3*NATOMS*IINTFP,Q)
      CALL PUTREC(20,'JOBARC','ORIENTMT',ININE*IINTFP,ORIENT)
      CALL PUTREC(20,'JOBARC','ATOMMASS',NATOMS*IINTFP,ATMASS)
      CALL PUTCREC(20, 'JOBARC', "PTGP    ", 4, FPGRP)
      CALL PUTCREC(20, 'JOBARC', "ABL_PTGP", 4, BPGRP)
      CALL PUTCREC(20, 'JOBARC', "CMP_PTGP", 4, PGRP) 
      call GETREC(0,'JOBARC','ZSYM',iTmp,indx)
C
      if (iTmp.lt.0) then
c      o write the atom symbols to JOBARC
c        (currently, this is only used by `xa2proc xyz`)
         indx = 1
         do i = 1, NATOMS
            ZSYMUNI(indx:indx+4) = ZSYM(i)(1:5)
            indx = indx + 5
         end do
         CALL PUTCREC(1,'JOBARC','ZSYM',5*NATOMS,ZSYMUNI)
      end if
      RETURN
      END
