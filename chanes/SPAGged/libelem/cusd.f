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
      SUBROUTINE CUSD1(IVAR)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*4 KLC,KLV,KNL,IVAR
      COMMON/MDLA/MT(15)

      MT(1)=2
      MT(2)=2
      MT(3)=2
      MT(4)=2
      KLC= 0
      KLV= 0
      KNL=4
      RETURN
      END


      SUBROUTINE CUSD2(OM,P1,L1,P2,L2,P3,L3)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/SUBC/Y(15,15),J(15)
      DOUBLE PRECISION P1    ,P2    ,P3    ,ZN
      DOUBLE COMPLEX Y,J
      DIMENSION P1(L1),P2(L2),P3(L3)

      DO 10 I=1,4
      DO 10 II=1,4
      Y(II,I)= DCMPLX(0.0D0,0.0D0)
  10  CONTINUE
      RETURN
      END


      SUBROUTINE CUSD3(NG,P1,L1,P2,L2,P3,L3,B1,KNC2,NR,*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION P1,P2,P3
      DIMENSION P1(L1),P2(L2),P3(L3)
      DOUBLE PRECISION B1
      DIMENSION B1(KNC2,NR)
      DOUBLE PRECISION LD
      UN1(U1,U2)=BT*(U1+VT)*(U1+VT)*(1.D0+LD*U2)*DTANH(AL*U2)
      BT=P3(1)
      LD=P3(2)
      AL=P3(3)
      VT=P3(4)

      DO 30 K=1,KNC2,2
      U1=B1(K,1)
      U2=B1(K,2)
      IF((U1+VT).LE.0.0D0)GOTO 25
      B1(K,1)=UN1(U1,U2)
      GOTO 20
  25  B1(K,1)=0.0D0
  20  B1(K+1,1)=0.0D0
  30  CONTINUE
      RETURN
      END


      SUBROUTINE CUSD4(NG,P1,L1,P2,L2,P3,L3,B1,KNC2,NR,*)

C
C

      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION P1,P2,P3
      DIMENSION P1(L1),P2(L2),P3(L3)
      DOUBLE PRECISION B1
      DIMENSION  B1(KNC2,NR)
      DOUBLE PRECISION LD
      TH(U)=DTANH(AL*U)
      DUN1(U1,U2)=2.D0*BT*(U1+VT)*(1.D0+LD*U2)*TH(U2)
      DUN2(U1,U2)=BT*(U1+VT)*(U1+VT)*((1.D0+LD*U2)*AL*(1.D0-TH(U2)*TH(U2
     +))+LD*TH(U2))
      BT=P3(1)
      LD=P3(2)
      AL=P3(3)
      VT=P3(4)
      DO 120 K=1,KNC2,2
      U1=B1(K,1)
      U2=B1(K,2)
      IF((U1+VT).LE.0.0D0) GOTO 125
      B1(K,1)=DUN1(U1,U2)
      B1(K,2)=DUN2(U1,U2)
      GOTO 130
 125  B1(K,1)=0.0D0
      B1(K,2)=0.0D0
 130  B1(K+1,1)=0.0D0
      B1(K+1,2)=0.0D0
 120  CONTINUE
      RETURN
      END


      SUBROUTINE CUSD5(NOI,NOU,EXIST,KOI,KOUV,KOPV,NR1V,NB1V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER NOI(5,2),NOU(5,2)
      LOGICAL EXIST(5,2)
      INTEGER KOI,KOUV(5),KOPV(5),NR1V(5),NB1V(5)
      NOI(1,1)=3
      NOI(1,2)=4
      NOU(1,1)=1
      NOU(1,2)=2
      NOU(2,1)=3
      NOU(2,2)=4
      EXIST(1,1)=.TRUE.
      EXIST(1,2)=.FALSE.
      EXIST(2,1)=.TRUE.
      EXIST(2,2)=.FALSE.
      KOI=1
      KOUV(1)=2
      KOPV(1)=0
      NR1V(1)=2
      NB1V(1)=2
      RETURN
      END


