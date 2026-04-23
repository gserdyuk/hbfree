!*==JUNC1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c





      SUBROUTINE junc1(Ivar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Mt
      INTEGER*4 klc , klv , knl , Ivar
      COMMON /mdla  / Mt(15)

      Mt(1) = 2
      Mt(2) = 2
      klc = 0
      klv = 0
      knl = 2
C     DEBUG SUBTRACE
      END


      SUBROUTINE junc2(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L1 , L2 , L3
      DOUBLE PRECISION Om
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Y , J
      DIMENSION P1(L1) , P2(L2) , P3(L3)

C  JUNC = DIODE WITHOUT NONLINEAR CAPACITANCE
      Y(1,1) = dcmplx(0.0D0,0.0D0)
      Y(1,2) = dcmplx(0.0D0,0.0D0)
      Y(2,1) = dcmplx(0.0D0,0.0D0)
      Y(2,2) = dcmplx(0.0D0,0.0D0)
C     DEBUG SUBTRACE
      END


      SUBROUTINE junc3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      ITEM MAT. MODELS OF SCHOTTKY CONTACT (WITHOUT CAPACITANCE)
C          PROBABLY HAVE TO INTRODUCE AN APPROXIMATION
C          AT HIGH CURRENTS. DON'T FORGET ABOUT JUNC4
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dun , u , um , un1
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION argmax/40.D0/
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION i0 , al
      dun(u) = i0*al*dexp(al*u)
      un1(u) = i0*(dexp(al*u)-1)
      i0 = P3(1)
      al = P3(2)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( al*u.LE.argmax ) THEN
            B1(k,1) = un1(u)
         ELSE
            um = argmax/al
            B1(k,1) = un1(um) + dun(um)*(u-um)
         ENDIF
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(P1,L1,P2,L2,P3,L3,KNC2,NR,I0,AL)
      END


      SUBROUTINE junc4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      SEE ITEM JUNC3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dun1 , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION i0 , al , argmin/ - 120.0D0/ , argmax/40.D0/
C                 $=MIN VALUE OF THE ARGUMENT OF THE EXP FUNCTION.
      dun1(u) = al*i0*dexp(al*u)
      i0 = P3(1)
      al = P3(2)

      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         B1(k,1) = 0.0D0
         IF ( u*al.GT.argmin .AND. u*al.LT.argmax ) B1(k,1) = dun1(u)
         IF ( u*al.GT.argmax ) B1(k,1) = dun1(argmax/al)
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE junc5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER kgr
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
      Noi(1,1) = 1
      Noi(1,2) = 2
      Nou(1,1) = 1
      Nou(1,2) = 2
      Exist(1,1) = .TRUE.
      Exist(1,2) = .FALSE.
      Koi = 1
      Kouv(1) = 1
      Kopv(1) = 0
      Nr1v(1) = 1
      Nb1v(1) = 1
      kgr = 1
      END


      SUBROUTINE junc6(Ng,P1,L1,P2,L2,P3,L3,Val,Dval,Kn,Nr,T)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Kn , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION T , zabs
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Val , Dval
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DIMENSION Val(Kn,Nr) , Dval(Kn,Nr)
      DOUBLE PRECISION uold , unew , ustep
      DOUBLE PRECISION delta , ubound , stepm

      delta = 0.01D0
      ubound = P3(3)
      stepm = P3(4)

      uold = 0.0D0
      unew = 0.0D0
      IF ( Kn.GE.2 ) THEN
         DO i = 2 , Kn
            uold = uold + zabs(Val(i,1))
            unew = unew + zabs(Val(i,1)+Dval(i,1))
         ENDDO
      ENDIF
C                            $ MAYBE NEED NOT + BUT -
      uold = uold + uold + dreal(Val(1,1))
      unew = unew + unew + dreal(Val(1,1)+Dval(1,1))
      ustep = unew - uold
C  MAYBE HERE $ NEED TO PUT ABS(     )

      IF ( unew.LE.ubound ) RETURN
      IF ( uold.LE.ubound ) THEN
         T = (ubound+delta-uold)/ustep
         RETURN
      ENDIF
      IF ( ustep.GT.stepm ) T = stepm/ustep
      RETURN
C     DEBUG SUBTRACE,INIT(UOLD,UNEW,UBOUND,STEPM,USTEP,T)
      END
