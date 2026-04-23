!*==MDSCH1.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c





      SUBROUTINE mdsch1(Ivar)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Mt
      INTEGER*4 klc , klv , knl , Ivar
      COMMON /mdla  / Mt(15)

      Mt(1) = Ivar
      Mt(2) = 2
      Mt(3) = 2
      klc = Ivar
      klv = 1 - klc
      knl = 2
C     DEBUG SUBTRACE
      END


      SUBROUTINE mdsch2(Om,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER L1 , L2 , L3
      DOUBLE PRECISION Om
      COMMON /subc  / Y(15,15) , J(15)
      DOUBLE PRECISION P1 , P2 , P3 , zn
      DOUBLE COMPLEX Y , J
      DIMENSION P1(L1) , P2(L2) , P3(L3)
C     WRITE(6,111) P3(1), P3(3), L3
C 111 FORMAT(2X,'MDSCH2: P3(1)=',E12.5,' P3(3)=',E12.5,' L3=',I3)
      zn = P3(1)**2 + (P3(3)*Om)**2
      Y(1,1) = dcmplx(P3(1)/zn,-Om*P3(3)/zn+P3(2)*Om)
      Y(1,2) = -dcmplx(0.0D0,Om*P3(2))
      Y(1,3) = -dcmplx(P3(1)/zn,-Om*P3(3)/zn)
      Y(2,1) = Y(1,2)
      Y(2,2) = dcmplx(0.0D0,Om*P3(2))
      Y(2,3) = dcmplx(0.0D0,0.0D0)
      Y(3,1) = Y(1,3)
      Y(3,2) = Y(2,3)
      Y(3,3) = dcmplx(P3(1)/zn,-Om*P3(3)/zn)
C     DEBUG SUBTRACE
      END


      SUBROUTINE mdsch3(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION eps2 , u , un1 , un2
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION fit , fr , tok , c0 , alfa , argmax/40.0D0/
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION B1
      un1(u) = tok*(dexp(alfa*u)-1.0D0)
      un2(u) = c0/dsqrt(1.0D0-u/fit)
      fit = P2(1)
      fr = P2(2)
      tok = P3(4)
      c0 = P3(5)
      alfa = P3(6)
      eps2 = 2.0D0*(fit-fr)
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*alfa.GT.argmax ) u = argmax/alfa
         IF ( u.GT.fr ) THEN
            B1(k,1) = un1(u) + un2(fr)*B1(k,2)*(1.0D0+(u-fr)/eps2)
         ELSE
            B1(k,1) = un1(u) + un2(u)*B1(k,2)
         ENDIF
         B1(k+1,1) = 0.0D0
      ENDDO
C     DEBUG SUBTRACE,INIT(FIT,FR,C0,TOK,ALFA,EPS2)
      END


      SUBROUTINE mdsch4(Ng,P1,L1,P2,L2,P3,L3,B1,Knc2,Nr,*)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DOUBLE PRECISION dun1 , fit2 , u , un2 , un3
      INTEGER k , Knc2 , L1 , L2 , L3 , Ng , Nr
      DOUBLE PRECISION P1 , P2 , P3
      DIMENSION P1(L1) , P2(L2) , P3(L3)
      DOUBLE PRECISION fit , fr , tok , c0 , alfa , argmax/40.0D0/
      DIMENSION B1(Knc2,Nr)
      DOUBLE PRECISION B1
      dun1(u) = tok*alfa*dexp(alfa*u)
      un2(u) = c0/dsqrt(1.0D0-u/fit)
      un3(u) = un2(u)/(2.0D0*(fit-u))

      fit = P2(1)
      fr = P2(2)
      tok = P3(4)
      c0 = P3(5)
      alfa = P3(6)
      fit2 = fit + fit
      DO k = 1 , Knc2 , 2
         u = B1(k,1)
         IF ( u*alfa.GT.argmax ) u = argmax/alfa
         IF ( u.GT.fr ) THEN
            B1(k,1) = dun1(u) + un3(fr)*B1(k,2)
            B1(k,2) = un2(fr) + un3(fr)*(u-fr)
         ELSE
            B1(k,1) = dun1(u) + un3(u)*B1(k,2)
            B1(k,2) = un2(u)
         ENDIF
         B1(k+1,1) = 0.0D0
         B1(k+1,2) = 0.0D0
C     IF(K.NE.KKK)GOTO 10
C     KKK=KKK+10
C     PRINT 111,(K,U,DUDT,B1(K,1),B1(K,2))
C 111 FORMAT('  K=',I4,' U=',E12.5,' DUDT=',E12.5,'B1(K,1)=',E12.5,
C    * ' B1(K,2)=',E12.5)
      ENDDO
C     DEBUG SUBTRACE
      END


      SUBROUTINE mdsch5(Noi,Nou,Exist,Koi,Kouv,Kopv,Nr1v,Nb1v)

      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER Noi(5,2) , Nou(5,2)
      LOGICAL Exist(5,2)
      INTEGER Koi , Kouv(5) , Kopv(5) , Nr1v(5) , Nb1v(5)

      Noi(1,1) = 3
      Noi(1,2) = 2
      Nou(1,1) = 3
      Nou(1,2) = 2
      Exist(1,1) = .TRUE.
      Exist(1,2) = .TRUE.
      Koi = 1
      Kouv(1) = 1
      Kopv(1) = 1
      Nr1v(1) = 2
      Nb1v(1) = 2

      END


      SUBROUTINE mdsch6(Ng,P1,L1,P2,L2,P3,L3,Val,Dval,Kn,Nr,T)
C
C
C
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
      ubound = P3(7)
      stepm = P3(8)

      uold = 0.0D0
      unew = 0.0D0
      IF ( Kn.GE.2 ) THEN
         DO i = 2 , Kn
            uold = uold + zabs(Val(i,1))
            unew = unew + zabs(Val(i,1)+Dval(i,1))
         ENDDO
      ENDIF

      uold = uold + uold + dreal(Val(1,1))
      unew = unew + unew + dreal(Val(1,1)+Dval(1,1))
      ustep = unew - uold
C  IT MAY BE NECESSARY TO USE ABS(     )

      IF ( unew.LE.ubound ) RETURN
      IF ( uold.LE.ubound ) THEN
         T = (ubound+delta-uold)/ustep
         RETURN
      ENDIF
      IF ( ustep.GT.stepm ) T = stepm/ustep
      RETURN
C     DEBUG SUBTRACE,INIT(USTEP,T)
      END
