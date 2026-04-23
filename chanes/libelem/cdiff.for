!*==CDIFF1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE cdiff1(Ivar)
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


      SUBROUTINE cdiff2(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L1 , L2 , L3
      DOUBLE PRECISION Om
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE PRECISION P1 , P2 , P3
      DOUBLE COMPLEX Y , J
      DIMENSION P1(L1) , P2(L2) , P3(L3)

      Y(1,1) = dcmplx(0.0D0,0.0D0)
      Y(1,2) = dcmplx(0.0D0,0.0D0)
      Y(2,1) = dcmplx(0.0D0,0.0D0)
      Y(2,2) = dcmplx(0.0D0,0.0D0)
C     DEBUG SUBTRACE
      END


      SUBROUTINE cdiff3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SEMICONDUCTOR MATERIAL MODEL OF DIFFUSION CAPACITANCE APPROXIMATION
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION tau , i0 , al , argmin/ - 120.0D0/ ,
     &                 argmax/40.0D0/

      c(u) = tau*al*i0*dexp(al*u)
      tau = P3(1)
      i0 = P3(2)
      al = P3(3)
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         B1(k,1) = 0.0D0
         IF ( u*al.GT.argmax ) u = argmax/al
         IF ( u*al.GT.argmin ) B1(k,1) = c(u)*(1+al*u)*B1(k,2)
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(P1,L1,P2,L2,P3,L3,KNC2,NR,I0,AL)
      END


      SUBROUTINE cdiff4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      SEE SEMICONDUCTOR CDIFF3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , tau , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION i0 , al , argmin/ - 120.0D0/ , argmax/40.0D0/
C                 $=MIN VALUE OF ARGUMENT OF EXP FUNCTION
      c(u) = tau*al*i0*dexp(al*u)
      tau = P3(1)
      i0 = P3(2)
      al = P3(3)
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*al.GT.argmax ) u = argmax/al
         IF ( u*al.LE.argmin ) THEN
            B1(k,1) = 0.D0
            B1(k,2) = 0.D0
         ELSE
            B1(k,1) = c(u)*al*(2.D0+al*u)*B1(k,2)
            B1(k,2) = c(u)*(1.D0+al*u)
         ENDIF
         B1(k+1,1) = 0.D0
         B1(k+1,2) = 0.D0
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE cdiff5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
      Noi(1,1) = 1
      Noi(1,2) = 2
      Nou(1,1) = 1
      Nou(1,2) = 2
      Exist(1,1) = .TRUE.
      Exist(1,2) = .TRUE.
      Koi = 1
      Kouv(1) = 1
      Kopv(1) = 1
      Nr1v(1) = 2
      Nb1v(1) = 2
      END


      SUBROUTINE cdiff6(Ng,P1,L1,P2,L2,P3,L3,Val,Dval,Kn,Nr,T)
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
C                            $ MAYBE NEED NOT + A -
      uold = uold + uold + dreal(Val(1,1))
      unew = unew + unew + dreal(Val(1,1)+Dval(1,1))
      ustep = unew - uold
C  MAYBE HERE $ NEED TO WRITE ABS(     )

      IF ( unew.LE.ubound ) RETURN
      IF ( uold.LE.ubound ) THEN
         T = (ubound+delta-uold)/ustep
         RETURN
      ENDIF
      IF ( ustep.GT.stepm ) T = stepm/ustep
      RETURN
C     DEBUG SUBTRACE,INIT(          USTEP,T)
      END
