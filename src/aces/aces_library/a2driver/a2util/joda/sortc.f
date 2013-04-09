
      SUBROUTINE SortC (NAtms, XX, AtNr, Y, NORD, X)
C     X is a scratch array of the same type and dimension as XX and Y
C     AtNr are the atomix numbers of the centers (dummy atom = 0)
C
C     SORTS VECTOR OF NUCLEAR COORDINATES - TO CHECK FOR EQUIVALENCE
C     OF TWO ORIENTATIONS - NEEDS Q VECTOR AND ATOMIC NUMBER VECTOR (AtNr)
C
      Implicit double precision (a-h,o-z)
      DIMENSION X(3*NATMS),XX(3*NATMS),Y(3*NATMS),NORD(NATMS)
      Integer AtNr(natms)
C
      ILINE(J)=1+J/3
C
C     SORT ON THE X - IF TWO X'S ARE EQUIVALENT, SORT ON Y AND SO ON.
C     
      DO 80 I=1,3*NATMS
 80      X(I)=XX(I)
C
C     FIRST GIVE DUMMY ATOMS RIDICULOUS SCRATCH COORDINATES - ENSURES
C     THAT THEY WILL WIND UP AT THE BOTTOM OF THE LIST
C
      DO 81 I=1,3*NATMS-2,3
         IF(AtNr(ILINE(I)) .eq. 0)THEN
            DO 82 J=0,2
 82            X(J+I) = -99995.
         ENDIF
 81   CONTINUE
      JK=1
 429  J=1
      DO 96 I=1,3*NATMS-2,3
C
C     CONTINUE WITH SORTING.
C
         IF(X(I)-X(J).GT.1D-6)J=I
         IF(DABS(X(I)-X(J)).LT.1D-6)THEN
            IF(X(I+1)-X(J+1).GT.1D-6)J=I
            IF(DABS(X(I+1)-X(J+1)).LT.1D-6)THEN
               IF(X(I+2)-X(J+2).GT.1D-6)J=I
            ENDIF
         ENDIF
 96   CONTINUE
      DO 93 I=0,2
C
C     Mass-WEIGHT SORTED VECTOR - WILL ZERO ELEMENTS CORRESPONDING
C     TO DUMMY ATOMS SO THEY DONT MUCK UP THE SYMMETRY CHECKING.
C     
         Y(3*JK-2+I)=X(J+I)*AtNr(ILINE(J))
 93      X(J)=-99999.D0
      NORD(JK)=(J+2)/3
      JK=JK+1
      if(jk.eq.NATMS+1)go to 999
      go to 429
 999  Continue
      Return
      end
