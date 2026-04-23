!*==CLIN1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE clin1(Ivar)
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


      SUBROUTINE clin2(Om,P1,L1,P2,L2,P3,L3)
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


      SUBROUTINE clin3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SEMICONDUCTOR MATERIAL MODEL OF BARRIER CAPACITANCE
C          =DEPENDENCE I(U)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION c0

      c(u) = c0
      c0 = P3(1)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         B1(k,1) = c(u)*B1(k,2)
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(P1,L1,P2,L2,P3,L3,KNC2,NR,I0,AL)
      END


      SUBROUTINE clin4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C      SEE SEMICONDUCTOR CLIN3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION c , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1 , c0
      DIMENSION B1(Knc2,Nr)

      c(u) = c0
      c0 = P3(1)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         B1(k,1) = 0.0D0
         B1(k,2) = c(u)
         B1(k+1,1) = 0.D0
         B1(k+1,2) = 0.D0
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE clin5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
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
