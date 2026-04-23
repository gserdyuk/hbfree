!*==CPOLY1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE cpoly1(Ivar)
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


      SUBROUTINE cpoly2(Om,P1,L1,P2,L2,P3,L3)
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


      SUBROUTINE cpoly3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SUBROUTINE FOR THE MATHEMATICAL MODEL OF BARRIER CAPACITANCE
C          =DEPENDENCY OF I(U)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , c1 , c2 , c3 , c4 , c5 , dcdu , u , u0
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION c0
      c(u) = c0 + ((((c5*u+c4)*u+c3)*u+c2)*u+c1)*u
      dcdu(u) = c1 + (((5.D0*c5*u+4.D0*c4)*u+3.D0*c3)*u+2.D0*c2)*u
      u0 = P3(1)
      c0 = P3(2)
      c1 = P3(3)
      c2 = P3(4)
      c3 = P3(5)
      c4 = P3(6)
      c5 = P3(7)
      DO k = 1 , Knc2 , 2
         u = B1(k,1) - u0
         B1(k,1) = (c(u)+dcdu(u)*u)*B1(k,2)
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(P1,L1,P2,L2,P3,L3,KNC2,NR,I0,AL)
      END


      SUBROUTINE cpoly4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SEE SUBROUTINE CPOLY3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , c1 , c2 , c3 , c4 , c5 , d2cdu2 , dcdu , u ,
     &                 u0
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1 , c0
      DIMENSION B1(Knc2,Nr)

      c(u) = c0 + ((((c5*u+c4)*u+c3)*u+c2)*u+c1)*u
      dcdu(u) = c1 + (((5.D0*c5*u+4.D0*c4)*u+3.D0*c3)*u+2.D0*c2)*u
      d2cdu2(u) = 2.D0*c2 + ((20.D0*c5*u+12.D0*c4)*u+6.D0*c3)*u
      u0 = P3(1)
      c0 = P3(2)
      c1 = P3(3)
      c2 = P3(4)
      c3 = P3(5)
      c4 = P3(6)
      c5 = P3(7)
      DO k = 1 , Knc2 , 2
         u = B1(k,1) - u0
         B1(k,1) = (dcdu(u)+dcdu(u)+u*d2cdu2(u))*B1(k,2)
         B1(k,2) = (c(u)+dcdu(u)*u)
         B1(k+1,1) = 0.D0
         B1(k+1,2) = 0.D0
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE cpoly5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
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
