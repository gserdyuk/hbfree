!*==CUSD1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



C************ ITUN DRAIN-SOURCE MODEL 1 PTBSH CURTICE ******************
C Quadratic approximation of the I-V characteristic for the field-effect transistor
C
C Isi(Uzi, Usi) = BETTA * (Uzi + VT)^2 * (1 + LAMDA * Usi) * TANH(ALF * Usi)
C
      SUBROUTINE cusd1(Ivar)
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
      END


      SUBROUTINE cusd2(Om,P1,L1,P2,L2,P3,L3)
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
      END


      SUBROUTINE cusd3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION al , bt , u1 , u2 , un1 , vt
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION ld
      un1(u1,u2) = bt*(u1+vt)*(u1+vt)*(1.D0+ld*u2)*dtanh(al*u2)
      bt = P3(1)
      ld = P3(2)
      al = P3(3)
      vt = P3(4)

      DO k = 1 , Knc2 , 2
         u1 = B1(k,1)
         u2 = B1(k,2)
         IF ( (u1+vt).LE.0.0D0 ) THEN
            B1(k,1) = 0.0D0
         ELSE
            B1(k,1) = un1(u1,u2)
         ENDIF
         B1(k+1,1) = 0.0D0
      ENDDO
      END


      SUBROUTINE cusd4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

C
C

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION al , bt , dun1 , dun2 , th , u , u1 , u2 , vt
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION ld
      th(u) = dtanh(al*u)
      dun1(u1,u2) = 2.D0*bt*(u1+vt)*(1.D0+ld*u2)*th(u2)
      dun2(u1,u2) = bt*(u1+vt)*(u1+vt)
     &              *((1.D0+ld*u2)*al*(1.D0-th(u2)*th(u2))+ld*th(u2))
      bt = P3(1)
      ld = P3(2)
      al = P3(3)
      vt = P3(4)
      DO k = 1 , Knc2 , 2
         u1 = B1(k,1)
         u2 = B1(k,2)
         IF ( (u1+vt).LE.0.0D0 ) THEN
            B1(k,1) = 0.0D0
            B1(k,2) = 0.0D0
         ELSE
            B1(k,1) = dun1(u1,u2)
            B1(k,2) = dun2(u1,u2)
         ENDIF
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
      ENDDO
      END


      SUBROUTINE cusd5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)
      Noi(1,1) = 3
      Noi(1,2) = 4
      Nou(1,1) = 1
      Nou(1,2) = 2
      Nou(2,1) = 3
      Nou(2,2) = 4
      Exist(1,1) = .TRUE.
      Exist(1,2) = .FALSE.
      Exist(2,1) = .TRUE.
      Exist(2,2) = .FALSE.
      Koi = 1
      Kouv(1) = 2
      Kopv(1) = 0
      Nr1v(1) = 2
      Nb1v(1) = 2
      END
