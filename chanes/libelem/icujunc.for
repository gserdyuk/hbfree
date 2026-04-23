!*==ICUJ1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




C******** ITUN ON VAX I=I0*ALPHA*(EXP(AL*U)-1.0) ******************
C    USED FOR MODELING ITUNs IN
C    BIPOLAR TRANSISTORS.
      SUBROUTINE icuj1(Ivar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Mt
      INTEGER*4 klc , klv , knl , Ivar
      COMMON /mdla  / Mt(15)

      Mt(1) = 2
      Mt(2) = 2
      Mt(3) = 2
      Mt(4) = 2
      klc = 0
      klv = 0
      knl = 4
C     DEBUG SUBTRACE
      END


      SUBROUTINE icuj2(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , ii , L1 , L2 , L3
      DOUBLE PRECISION Om
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Y , J
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      DO i = 1 , 4
         DO ii = 1 , 4
            Y(ii,i) = dcmplx(0.0D0,0.0D0)
         ENDDO
      ENDDO

C     DEBUG SUBTRACE
      END


      SUBROUTINE icuj3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SUBROUTINE FOR THE MATHEMATICAL MODEL OF A BIPOLAR TRANSISTOR ITUN.
C          PROBABLY NEEDS AN APPROXIMATION
C          AT HIGH CURRENTS.
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION u , un1
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION i0 , al , alpha , argmax/40.0D0/

      un1(u) = alpha*i0*(dexp(al*u)-1)
      alpha = P3(1)
      i0 = P3(2)
      al = P3(3)

      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( al*u.GT.argmax ) u = argmax/al
         B1(k,1) = un1(u)
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(P1,L1,P2,L2,P3,L3,KNC2,NR,I0,AL,ALPHA,U)
      END


      SUBROUTINE icuj4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      SEE SUBROUTINE ICUJ3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dun1 , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION i0 , al , alpha , argmin/ - 120.0D0/ ,
     &                 argmax/40.0D0/
C                      $=MIN VALUE OF THE ARGUMENT OF THE EXP FUNCTION.

      dun1(u) = alpha*al*i0*dexp(al*u)
      alpha = P3(1)
      i0 = P3(2)
      al = P3(3)

      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( al*u.GT.argmax ) u = argmax/al
         B1(k,1) = 0.0D0
         IF ( u*al.GT.argmin ) B1(k,1) = dun1(u)
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE icuj5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER kgr
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
      Noi(1,1) = 3
      Noi(1,2) = 4
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


      SUBROUTINE icuj6(Ng,P1,L1,P2,L2,P3,L3,Val,Dval,Kn,Nr,T)
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
      ubound = P3(4)
      stepm = P3(5)

      uold = 0.0D0
      unew = 0.0D0
      IF ( Kn.GE.2 ) THEN
         DO i = 2 , Kn
            uold = uold + zabs(Val(i,1))
            unew = unew + zabs(Val(i,1)+Dval(i,1))
         ENDDO
      ENDIF
C                            $ MAY NEED NOT + BUT -
      uold = uold + uold + dreal(Val(1,1))
      unew = unew + unew + dreal(Val(1,1)+Dval(1,1))
      ustep = unew - uold
C  MAYBE HERE $ NEEDS TO BE REPLACED WITH ABS(     )

      IF ( unew.LE.ubound ) RETURN
      IF ( uold.LE.ubound ) THEN
         T = (ubound+delta-uold)/ustep
         RETURN
      ENDIF
      IF ( ustep.GT.stepm ) T = stepm/ustep
      RETURN
C     DEBUG SUBTRACE,INIT(          USTEP,T)
      END
