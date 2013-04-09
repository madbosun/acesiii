      SUBROUTINE GEN_NEWXYZFRMOLDQ
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
#include "mxatms.par"
#include "coord.com"
C
      LOGICAL OPTRES
      INTEGER NA(MXATMS),NB(MXATMS),NC(MXATMS)
      DIMENSION A(3,MXATMS),Q0(3*MXATMS)
C
      COMMON /USINT/ NX, NXM6, IARCH, NCYCLE, NUNIQUE, NOPT
      COMMON /FLAGS/ IFLAGS(100),IFLAGS2(500)
      COMMON /OPTCTL/ IPRNT,INR,IVEC,IDIE,ICURVY,IMXSTP,ISTCRT,IVIB,
     $   ICONTL,IRECAL,INTTYP,IDISFD,IGRDFD,ICNTYP,ISYM,IBASIS,
     $   XYZTol
      COMMON/RESTART2/OPTRES
C
      LOGICAL LINEAR
#include "io_units.par"
C
      NUMAT=0
      CCOS=0.D0
      CALL ZERO(Q,3*MXATMS)
      IONE=1
      CALL GETREC(-1,'JOBARC','PASS1   ',IONE,IX)
      NDIS=0
      If (IPrnt .ge. 20) then
         Write (LuOut, *) 'GMETRY starting with R'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf
      IF ( ncycle .NE. 0. or. ix.ne.0 ) THEN
C
C     DECOMPRESS R
C
 
         If (IPrnt .ge. 20) Write (LuOut, '(a)')
     $      '@GMETRY-I, Decompressing R.'
         CALL USQUSH(R,NXM6)
         IF(.NOT.OPTRES)  WRITE(LuOut,78)
 78      FORMAT(' Updating structure...')
      ELSE
         ATOR=DACOS(-1.D0)/180.D0
         ATOB=0.529177249D0
         If (IPrnt .ge. 20) Write (LuOut,*)
     $      'GMETRY: Converting to radians & bohr.',ator,atob
         DO 19 IX=8,NX-1,3
            IF(IX.NE.8)R(IX+1)=R(IX+1)*ATOR
 19         R(IX)=R(IX)*ATOR
         IF(IFLAGS(78).EQ.0)THEN
          DO 18 IX=4,NX-2,3
 18          R(IX)=R(IX)/ATOB
         ENDIF
      ENDIF
C
      If (IPrnt .ge. 20) then
         Write (LuOut, *) 'GMETRY using R vector'
         Write (LuOut,'((I3,1X,f12.6))') (i,R(i),i=1,NX)
      EndIf
C
      DO 1 I=1,NATOMS
         DO 2 J=1,3
            A(J,I)=R(3*I+J-3)
 2       CONTINUE
         NA(I)=NCON(3*I)
         NB(I)=NCON(3*I-1)
         NC(I)=NCON(3*I-2)
 1    CONTINUE
      LINEAR=.FALSE.
C
C     *** FIND THE FIRST ATOM (IF ANY) WHICH IS NOT COLLINEAR WITH ATOMS
C     1 AND 2.  ASSUME THAT THE ATOMS ARE COLLINEAR IF THE COSINE OF THE
C     ANGLE IS GREATER THAN 0.99999
C
      DO 10 I=3,NATOMS
         CCOS=DCOS(A(2,I))
         IF (DABS(CCOS).LT.0.999999) GO TO 20
 10   CONTINUE
C
C     If we get here, the molecule is linear
C
      IF(NATOMS.GT.2)THEN
       WRITE(LUOUT,8080)
       WRITE(LUOUT,8081)
8080   FORMAT(T3,'@GMETRY-F, Linear frameworks not permitted.')
8081   FORMAT(T3,'           Try using dummy atoms in the Z-matrix.')
       CALL ERREX
      ENDIF
      If (IPrnt .ge. 20) Write (LuOut,*) 'Molecule is linear.'
      LINEAR=.TRUE.
      I=NATOMS+1
C
 20   Q(3*I-5)=A(1,I-1)
      Q(3*I-4)=0.0D0
      Q(3*I-3)=0.0D0
      II=I-2
      DO 30 J=1,II
         K=II-J+1
         Q(3*K-2)=Q(3*K+1)-A(1,K+1)
         Q(3*K-1)=0.0D0
         Q(3*K)=0.0D0
 30   CONTINUE
C
      IF (.NOT. LINEAR) THEN
         SSIN=DSIN(A(2,I))
         Q(3*I-2)=Q(3*I-5)-A(1,I)*CCOS
         Q(3*I-1)=A(1,I)*SSIN
         Q(3*I)=0.0D0
         IF (I.EQ.NATOMS) GO TO 80
         II=I+1
         DO 70 I=II,NATOMS
            COSA=DCOS(A(2,I))
            MB=NB(I)
            MC=NC(I)
            XB=Q(3*MB-2)-Q(3*MC-2)
            YB=Q(3*MB-1)-Q(3*MC-1)
            ZB=Q(3*MB)-Q(3*MC)
            RBC=1.0D0/DSQRT(XB*XB+YB*YB+ZB*ZB)
            IF (ABS(COSA) .GE. 0.999999) THEN
C
C     ATOMS MC, MB, AND (I) ARE COLLINEAR
C
               RBC=A(1,I)*RBC
               JX=3*I-2
               Q(JX)=Q(3*MC-2)-XB*RBC
               Q(JX+1)=Q(3*MC-1)-YB*RBC
               Q(JX+2)=Q(3*MC)-ZB*RBC
            ELSE
C
C     THE ATOMS ARE NOT COLLINEAR
C
               MA=NA(I)
               XA=Q(3*MA-2)-Q(3*MC-2)
               YA=Q(3*MA-1)-Q(3*MC-1)
               ZA=Q(3*MA)-Q(3*MC)
C
C     ROTATE ABOUT THE Z-AXIS TO MAKE YB=0, AND XB POSITIVE.  IF XYB IS
C     TOO SMALL, FIRST ROTATE THE Y-AXIS BY 90 DEGREES.
C
               XYB=DSQRT(XB*XB+YB*YB)
               K=-1
               IF (XYB .LE. 0.1) THEN
                  XPA=ZA
                  ZA=-XA
                  XA=XPA
                  XPB=ZB
                  ZB=-XB
                  XB=XPB
                  XYB=DSQRT(XB*XB+YB*YB)
                  K=+1
               ENDIF
C
C     ROTATE ABOUT THE Y-AXIS TO MAKE ZB VANISH
C
               COSTH=XB/XYB
               SINTH=YB/XYB
               XPA=XA*COSTH+YA*SINTH
               YPA=YA*COSTH-XA*SINTH
               SINPH=ZB*RBC
               COSPH=DSQRT(DABS(1.0-SINPH*SINPH))
               XQA=XPA*COSPH+ZA*SINPH
               ZQA=ZA*COSPH-XPA*SINPH
C
C     ROTATE ABOUT THE X-AXIS TO MAKE ZA=0, AND YA POSITIVE.
C
               YZA=DSQRT(YPA**2+ZQA**2)
               COSKH=YPA/YZA
               SINKH=ZQA/YZA
C
C     COORDINATES :-   A=(XQA,YZA,0),   B=(RBC,0,0),  C=(0,0,0)
C     NONE ARE NEGATIVE.
C     THE COORDINATES OF I ARE EVALUATED IN THE NEW FRAME.
C
               SINA=DSIN(A(2,I))
               SIND=DSIN(A(3,I))
               COSD=DCOS(A(3,I))
               XD=A(1,I)*COSA
               YD=A(1,I)*SINA*COSD
               ZD=A(1,I)*SINA*SIND
C
C     TRANSFORM THE COORDINATES BACK TO THE ORIGIONAL SYSTEM.
C
               YPD=YD*COSKH-ZD*SINKH
               ZPD=ZD*COSKH+YD*SINKH
               XPD=XD*COSPH-ZPD*SINPH
               ZQD=ZPD*COSPH+XD*SINPH
               XQD=XPD*COSTH-YPD*SINTH
               YQD=YPD*COSTH+XPD*SINTH
               IF (K .GE. 1) THEN
                  XRD=-ZQD
                  ZQD=XQD
                  XQD=XRD
               ENDIF
 60            Q(3*I-2)=XQD+Q(3*MC-2)
               Q(3*I-1)=YQD+Q(3*MC-1)
               Q(3*I)=ZQD+Q(3*MC)
            ENDIF
 70      CONTINUE
C
C     *** NOW REMOVE THE DUMMY ATOM COORDINATES, IF ANY, FROM THE ARRAY
C
C
C     PRDDO REQUIRES INPUT AS ATOMIC UNITS
C
      ENDIF
 80   SCALE=1.0D0
      If (IPrnt .ge. 200) then
        Write (LuOut, 444)
 444    Format(T3,'@GMETRY-I, Cartesian coordinates before scaling:')
         Write (LuOut, 445) (i,Q(3*i-2),Q(3*i-1),Q(3*i),
     $      i=1,NAtoms)
 445    Format(T5,I3,1X,F12.6,1X,F12.6,1X,F12.6)
      EndIf
      IF (NDIS.NE.0) SCALE=0.529177249D0
      DO 115 I=1,NATOMS
         Q(3*I-2)=Q(3*I-2)/SCALE
         Q(3*I-1)=Q(3*I-1)/SCALE
 115     Q(3*I)=Q(3*I)/SCALE
      If (IPrnt .ge. 200) then
        Write (LuOut, 446)
 446    Format(T3,'@GMETRY-I, Cartesian coordinates after scaling:')
         Write (LuOut, 445) (i,Q(3*i-2),Q(3*i-1),Q(3*i),
     $      i=1,NAtoms)
      EndIf
      IF (NATOMS.EQ.NUMAT) RETURN
      J=0
      DO 100 I=1,NATOMS
         IF (NR(I) .NE. 99) THEN
            J=J+1
C
C     CHANGE TO CONVENTION IN GAUSSIAN.  FOR DEBUGGING PURPOSES ONLY.
C
C           Q(3*J-2)=Q(3*I-2)
C           Q(3*J-1)=Q(3*I-1)
C           Q(3*J)=Q(3*I)
            Q0(3*J-2)=Q(3*I-2)
            Q0(3*J-1)=Q(3*I-1)
            Q0(3*J)=Q(3*I)
            Q(3*J)=Q0(3*I-2)
            Q(3*J-2)=Q0(3*I-1)
            Q(3*J-1)=-Q0(3*I)
         ENDIF
 100  CONTINUE
      If (IPrnt .ge. 20) then
         Write (LuOut, 555)
 555     Format(T3,'@GMETRY-I, Cartesian coordinates from Z-matrix:')
         Write (LuOut, 445) (i,Q(3*i-2),Q(3*i-1),Q(3*i),
     $      i=1,NAtoms)
      EndIf
C C
C DEAL WITH CASE WHERE ANGLES OR DIHEDRALS NOT WITHIN BOUNDS HERE.
C
      CALL ZERO(Q0,3*MXATMS)
      CALL VADD(Q0,R,Q0,NX,1.D0)
      CALL XTOR(R,0)
      CALL USQUSH(R,NXM6)
      DO 107 I=1,NX
       IF(DABS(DABS(R(I))-DABS(Q0(I))).GT.1.D-4)THEN
        WRITE(6,8331)I,Q0(I),R(I)
 8331   FORMAT('@GMETRY-W, Internal coordinate #',i3,
     &' not within bounds:',
     &/,t3,' Value was: ',f10.5,' and has been changed to ',f10.5,'.')
       ENDIF
 107  CONTINUE
      RETURN
      END
