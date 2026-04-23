!*==CBARR1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE cbarr1(Ivar)
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


      SUBROUTINE cbarr2(Om,P1,L1,P2,L2,P3,L3)
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


      SUBROUTINE cbarr3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SEMICONDUCTOR MATERIAL MODEL OF BARRIER CAPACITANCE
C          =DEPENDENCE I(U)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , f1 , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION c0 , fi0 , an , ug

      f1(u) = an/(fi0-u)
      c(u) = c0/(1.D0-u/fi0)**an
      c0 = P3(1)
      fi0 = P3(2)
      an = P3(3)
      ug = P3(4)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u.LT.ug ) THEN
            B1(k,1) = c(u)*(1.D0+f1(u)*u)*B1(k,2)
         ELSE
            B1(k,1) = c(ug)*(1.D0+f1(ug)*(2.D0*u-ug))*B1(k,2)
         ENDIF
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(C0,FI0,AN,UG,U,B1  )
      END


      SUBROUTINE cbarr4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      SEE SEMICONDUCTOR CBARR3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , f1 , f2 , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1 , c0 , fi0 , an , ug
      DIMENSION B1(Knc2,Nr)
      f1(u) = an/(fi0-u)
      c(u) = c0/(1-u/fi0)**an
      f2(u) = (an+1.D0)/(fi0-u)
      c0 = P3(1)
      fi0 = P3(2)
      an = P3(3)
      ug = P3(4)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u.LT.ug ) THEN
            B1(k,1) = c(u)*f1(u)*(2.D0+f2(u)*u)*B1(k,2)
            B1(k,2) = c(u)*(1.D0+f1(u)*u)
            B1(k+1,1) = 0.D0
            B1(k+1,2) = 0.D0
         ELSE
            B1(k,1) = 2.D0*c(ug)*f1(ug)*B1(k,2)
            B1(k,2) = c(ug)*(1.D0+f1(ug)*(2.D0*u-ug))
            B1(k+1,1) = 0.D0
            B1(k+1,2) = 0.D0
         ENDIF
      ENDDO
C     DEBUG INIT(C0,AN,FI0,UG,B1,U  )
      END


      SUBROUTINE cbarr5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
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
C
      END
