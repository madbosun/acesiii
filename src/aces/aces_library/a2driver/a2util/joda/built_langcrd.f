      SUBROUTINE BULT_LANGCRD(CARTCOORD, BMATRX, ANGL, ICON1, ICON2, 
     &                        ICON3, ICON4, IANGS, NRATMS, TOTREDNCO,
     &                        EPSILON)
C
C Setup the angle bending B-matrix elements for linear bonds. In this case
C we must define two angle bend coordinates.
C
C Assuming that the projection is done on XY-plane. The Z-axis is the 
C perpendicular axis. Given r_i = r_ba and r_j = r_bc then SinG_ij =  
C (x_1*y_2 - y_1*x_2)/R_i*R_j
C
C B(*,a,*,*)(x) =0.0, B(*,*,*,c)(x) = 0.0, 
C B(*,*,b,*) = B(*,a,*,*) - B(*,*,*,c)
C B(*,a,*,*)(y) = -(z_a-z_b)/R_BA, B(*,a,*,*)(z) =  (y_a-y_b)/R_BA
C B(*,*,*,c)(y) =  (z_c-z_b)/R_BC, B(*,*,*,c)(z) = -(y_c-y_b)/R_BC
C B(*,*,b,*)(y) = B(*,a,*,*)(y) + B(*,*,*,c)(y)
C B(*,*,b,*)(z) = B(*,a,*,*)(z) + B(*,*,*,c)(z)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      INTEGER TOTREDNCO
      DIMENSION CARTCOORD(3*NRATMS), BMATRX(TOTREDNCO, 3*NRATMS),
     &          VECBA(3), VECBC(3), VECAC(3), VECX(3), VECY(3), 
     &          VECZ(3), VECBCX(3), VECBAX(3),  VECN2AC(3),
     &          VECBN2A(3), VECBN2B(3), VECW(3)

C
      DATA VECX /1,0,0/, VECY /0,1,0/, VECZ /0,0,1/
      DATA ZERO /0.0D0/
      DATA MONE /-1/
C

          Write(6,*)
          Write(6, "(a)") "Entering the special linear angle block"

C
      DISTBA = DIST(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON2 - 2))
      DISTBC = DIST(CARTCOORD(3*ICON3 - 2), CARTCOORD(3*ICON2 - 2))

      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON1 - 2), 
     &         VECBA, 1)
      CALL VEC(CARTCOORD(3*ICON2 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECBC, 1)
      CALL VEC(CARTCOORD(3*ICON1 - 2), CARTCOORD(3*ICON3 - 2), 
     &         VECAC, 1)     
C
      CALL NORMAL(VECBA, 3)
      CALL NORMAL(VECBC, 3)
      CALL CROSS(VECBA, VECBC, VECN2AC, 0)
      CALL NORMAL(VECN2AC, 3)


          Write(6, *)
          Write(6, *) "The normal vector"
          Write(6, "(a, 3F10.4)") "VECN2AC", (VECN2AC(I), I=1,3)

C
      CALL CROSS(VECBC, VECN2AC, VECW, 1)
      CALL CROSS(VECBC, VECW, VECBN2A, 1)
      CALL CROSS(VECW, VECBA, VECBN2B, 1)
C

          Write(6, *)
          Write(6, *) "The basis vctors for special angles"
          Write(6, "(a, 3F10.4)") "VECBN2A", (VECBN2A(I), I=1,3)
          Write(6, "(a, 3F10.4)") "VECBN2B", (VECBN2B(I), I=1,3)
          Write(6,*)

C
      DO IXYZ = 1, 3
         BMATRX(IANGS, (3*ICON3 - 3) + IXYZ) =
     &                                 VECBN2A(IXYZ)/DISTBC
         BMATRX(IANGS, (3*ICON2 - 3) + IXYZ) =
     &                               - VECBN2A(IXYZ)/DISTBC
     &                               - VECBN2B(IXYZ)/DISTBA
         BMATRX(IANGS, (3*ICON1 - 3) + IXYZ) =
     &                                 VECBN2B(IXYZ)/DISTBA
      ENDDO
C
C Assign the perpendicular bond angle
C
      DINVPI = (ATAN(DFLOAT(1))*DFLOAT(4))/180.0D0
      ANGL = 180.0D0*DINVPI
C
C----------This block was old an old algorithm----------
C----------The old block end here ----------
C
      RETURN
      END
         
