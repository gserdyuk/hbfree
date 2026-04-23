!*==ICUPL1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



C******** ITUN ON VAX I=I0+....+I5(U-U0)**5 ******************
C    USED FOR MODELING ITUN'S
C    ANYWHERE.
      SUBROUTINE icupl1(Ivar)
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


      SUBROUTINE icupl2(Om,P1,L1,P2,L2,P3,L3)
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


      SUBROUTINE icupl3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SUBROUTINE FOR THE MATHEMATICAL MODEL
C      OF A POLYNOMIAL ITUN.
C
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION u , un1
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION a0 , a1 , a2 , a3 , a4 , a5 , u0
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION B1
      un1(u) = a0 + u*(a1+u*(a2+u*(a3+u*(a4+u*a5))))
      u0 = P3(1)
      a0 = P3(2)
      a1 = P3(3)
      a2 = P3(4)
      a3 = P3(5)
      a4 = P3(6)
      a5 = P3(7)
C
      DO k = 1 , Knc2 , 2
         u = B1(k,1) - u0
         B1(k,1) = un1(u)
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT
      END


      SUBROUTINE icupl4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
C      SEE SUBROUTINE ICUJ3
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dun1 , u
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION B1
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION a1 , a2 , a3 , a4 , a5 , u0
      DIMENSION B1(Knc2,Nr)
      dun1(u) = a1 + u*(2.D0*a2+u*(3.D0*a3+u*(4.D0*a4+u*5.D0*a5)))
      u0 = P3(1)
      a1 = P3(3)
      a2 = P3(4)
      a3 = P3(5)
      a4 = P3(6)
      a5 = P3(7)

      DO k = 1 , Knc2 , 2
         u = B1(k,1) - u0
         B1(k,1) = dun1(u)
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE
      END

      SUBROUTINE icupl5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
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
