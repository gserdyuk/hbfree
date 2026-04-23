c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c



      SUBROUTINE    CURT1(IVAR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER*4     MT
      INTEGER*4     KLC,KLV,KNL
C
      COMMON/MDLA/  MT(15)
C
      MT(1)=IVAR
      MT(2)=IVAR
      MT(3)=IVAR
      MT(4)=2
      MT(5)=2
      MT(6)=2
C
      KLC=(1-IVAR)*4
      KLV=4-KLC
      KNL=4
C
      RETURN
      END
      SUBROUTINE    CURT2(OM,P1,L1,P2,L2,P3,L3)
C
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION          P1,P2,P3
      DOUBLE COMPLEX       Y,J,ZERO/(0.0D0,0.0D0)/,CR,CI
C
      COMMON/SUBC/  Y(15,15),J(15)
C
      DIMENSION     P1(L1),P2(L2),P3(L3)
C
      CR(X)=DCMPLX(X,0.0D0)
      CI(X)=DCMPLX(0.0D0,OM*X)
C
      G1=1.D0/P3(1)
      G2=1.D0/P3(2)
      G3=1.D0/P3(3)
      IF(OM.EQ.0.D0) THEN
      G31=0.000001D0
      G23=0.000001D0
      ELSE
      G31=1.D0/P3(14)
      G23=1.D0/P3(15)
      ENDIF
      C12=P3(4)
      C13=P3(5)
C
      DO 10 I=1,8
      DO 10 K=1,8
   10 Y(K,I)=ZERO
C
      Y(1,1)=CR(G2)
      Y(1,5)=-CR(G2)
      Y(2,2)=CR(G1)
      Y(2,4)=-CR(G1)
      Y(3,3)=CR(G3)
      Y(3,6)=-CR(G3)
      Y(4,2)= Y(2,4)
      Y(4,4)=CR(G1)+CI(C12)+CI(C13)+CR(G31)
      Y(4,5)=-CI(C12)
      Y(4,6)=-CI(C13)-CR(G31)
      Y(5,1)= Y(1,5)
      Y(5,4)= Y(4,5)
      Y(5,5)=CR(G2)+CI(C12)+CR(G23)
      Y(5,6)=-CR(G23)
      Y(6,5)= Y(5,6)
      Y(6,3)= Y(3,6)
      Y(6,4)= Y(4,6)
      Y(6,6)= CR(G3)+CI(C13)

C
      RETURN
      END
      SUBROUTINE      CURT3(NG,P1,L1,P2,L2,P3,L3,                      B
     +1,KNC2,NR,*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION      B1(KNC2,NR)
       DIMENSION      P1(L1),P2(L2),P3(L3)
C
       DOUBLE PRECISION           P1,P2,P3
       DOUBLE PRECISION           IS,LD
       DOUBLE PRECISION           B1
       DOUBLE PRECISION UEPS/0.05D0/,ARGMAX/40.0D0/
C
      TH(U)=DTANH(AL*U)
      UN1(U)=IS*(DEXP(AS*U)-1.D0)
      UN2(U)=C230/DSQRT(1.0D0-U/VBI)
      UN3(U1,U2)=BT*(U1+VT)*(U1+VT)*(1.D0+LD*U2)*TH(U2)
C
      AL=P3(6)
      BT=P3(7)
      LD=P3(8)
      VT=P3(9)
      IS=P3(10)
      AS=P3(11)
      C230=P3(12)
      VBI=P3(13)
C     UEPS=0.05
      VBI0=VBI-UEPS
C
      IF(NG.NE.1) GOTO 100
C
      DO 10 K=1,KNC2,2
      U=B1(K,1)
      IF(U*AS.GT.ARGMAX) U=ARGMAX/AS
      IF(U.GT.VBI0) GOTO 5
C     THE DISPLACEMENT CURRENT IS DEFINED AS C(U)*DU/DT AND NOT
C     D(C(U)*U)/DT=(C(U)+U*DC/DU)*DU/DT
C
      B1(K,1)=UN1(U)+UN2(U)*B1(K,2)
      GOTO 9
C
C     IF U > VBI0 = VBI - UEPS, THEN THE CAPACITANCE IS CALCULATED BY
C     THE FORMULA C(U) = C(VBI0) + DC/DU(U=VBI0) * (U - VBI0)
C
    5 B1(K,1)=UN1(U)+UN2(VBI0)*(1.D0+(U-VBI0)/(2.D0*UEPS))*B1(K,2)
C
    9 B1(K+1,1)=0.0D0

   10 CONTINUE
      RETURN
C
  100 CONTINUE
      DO 30 K=1,KNC2,2
      U1=B1(K,1)
      U2=B1(K,2)
C
      IF((U1+VT).LE.0.0D0)GOTO 25
      B1(K,1)=UN3(U1,U2)
      GOTO 20
   25 B1(K,1)=0.0D0
   20 B1(K+1,1)=0.0D0
   30 CONTINUE
      RETURN
C
      END

      SUBROUTINE      CURT4(NG,P1,L1,P2,L2,P3,L3,                      B
     +1,KNC2,NR,*)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION      P1(L1),P2(L2),P3(L3)
      DIMENSION      B1(KNC2,NR)
C
      DOUBLE PRECISION           P1,P2,P3
      DOUBLE PRECISION           B1
      DOUBLE PRECISION           IS,LD
      DOUBLE PRECISION   UEPS/0.05D0/,ARGMIN/-120.0D0/,ARGMAX/40.0D0/
      DOUBLE PRECISION   SMIN/1.0D-04/,GMIN/1.0D-04/
C
      TH(U)=DTANH(AL*U)
      DUN1(U)=IS*AS*DEXP(AS*U)
      UN2(U)=C230/DSQRT(1.0D0-U/VBI)
      UN3(U)=UN2(U)/2*(VBI-U)
      DUN2(U1,U2)=2.D0*BT*(U1+VT)*(1.D0+LD*U2)*TH(U2)
      DUN3(U1,U2)=BT*(U1+VT)*(U1+VT)*((1.D0+LD*U2)*AL*(1.D0-TH(U2)*TH(U2
     +))+LD*TH(U2))
C
      AL=P3(6)
      BT=P3(7)
      LD=P3(8)
      VT=P3(9)
      IS=P3(10)
      AS=P3(11)
      C230=P3(12)
      VBI=P3(13)
C     UEPS=0.05
      VBI0=VBI-UEPS
C
      IF(NG.NE.1) GOTO 100
      DO 10 K=1,KNC2,2
      U=B1(K,1)
      IF(U*AS.GT.ARGMAX) U=ARGMAX/AS
      IF(U.GT.VBI0) GOTO 3
      IF(U*AS.LT.ARGMIN) GO TO 2
C  NORMAL CASE:
      B1(K,1)=DUN1(U)+UN3(U)*B1(K,2)
      B1(K,2)=UN2(U)
      GOTO 5
C  U IS TOO SMALL ( EXP(AS*U) <= 1.E-60 )
    2 B1(K,1)=UN3(U)*B1(K,2)
      B1(K,2)=UN2(U)
      GOTO 5
C  U IS CLOSE TO VBI. UN2(U) -> K SINGULARITY
    3 CONTINUE
      B1(K,1)=DUN1(U)+UN2(VBI0)/(2.D0*UEPS)*B1(K,2)
      B1(K,2)=UN2(VBI0)*(1.D0+(U-VBI0)/(2.D0*UEPS))
C
    5 CONTINUE
      B1(K+1,1)=0.0D0
      B1(K+1,2)=0.0D0
C
   10 CONTINUE
      RETURN
C
  100 CONTINUE
      DO 20 K=1,KNC2,2
      U1=B1(K,1)
      U2=B1(K,2)
      IF((U1+VT).LE.0.0D0) GOTO 25
      B1(K,1)=DUN2(U1,U2)
      B1(K,2)=DUN3(U1,U2)
      GOTO 30
C
   25 CONTINUE
      B1(K,1)=0.0D0
      B1(K,2)=0.0D0
C
   30 CONTINUE
      B1(K+1,1)=0.0D0
      B1(K+1,2)=0.0D0
   20 CONTINUE
C
      RETURN
C
      END
      SUBROUTINE     CURT5(NOI,NOU,EXIST,KOI,KOUV,                     K
     +OPV,NR1V,NB1V)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       LOGICAL       EXIST(5,2)
       INTEGER       NOI(5,2),NOU(5,2)
       INTEGER       KOI,KOUV(5),KOPV(5),NR1V(5),NB1V(5)
C
      KOI=2
C
      NOI(1,1)=5
      NOI(1,2)=6
C
      NOU(1,1)=5
      NOU(1,2)=6
C
      EXIST(1,1)=.TRUE.
      EXIST(1,2)=.TRUE.
C
      KOUV(1)=1
      KOPV(1)=1
C
      NR1V(1)=2
      NB1V(1)=2
C
      NOI(2,1)=4
      NOI(2,2)=6
C
      NOU(2,1)=5
      NOU(2,2)=6
      NOU(3,1)=4
      NOU(3,2)=6
C
      KOUV(2)=2
      KOPV(2)=0
C
      NR1V(2)=2
      NB1V(2)=2
C
      EXIST(2,1)=.TRUE.
      EXIST(2,2)=.FALSE.
      EXIST(3,1)=.TRUE.
      EXIST(3,2)=.FALSE.
      RETURN
C
      END
      SUBROUTINE CURT6(NG,P1,L1,P2,L2,P3,L3,VAL,DVAL,                 KN
     +,NR,T)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION P1,P2,P3
      DOUBLE COMPLEX VAL,DVAL
      DIMENSION P1(L1),P2(L2),P3(L3)
      DIMENSION VAL(KN,NR),DVAL(KN,NR)
      DOUBLE PRECISION UOLD,UNEW,USTEP
      DOUBLE PRECISION DELTA,UBOUND,STEPM

      DELTA=0.01D0
C      UBOUND=0.6
C      STEPM =0.3
       UBOUND=P3(16)
       STEPM=P3(17)
C      UBOUND י STEPM גשלי 0.5 י 0.1 .
      IF(NG.NE.1)GOTO 100
      UOLD=0.0D0
      UNEW=0.0D0
      IF(KN.LT.2)GO TO 15
      DO 10 I=2,KN
      UOLD=UOLD+ZABS(VAL(I,1))
   10 UNEW=UNEW+ZABS(VAL(I,1)+DVAL(I,1))
   15 CONTINUE

      UOLD=UOLD+UOLD+DREAL(VAL(1,1))
      UOLD1=UOLD
      UNEW=UNEW+UNEW+DREAL(VAL(1,1)+DVAL(1,1))
      UNEW1=UNEW
      USTEP=    UNEW-UOLD
      IF(UNEW.LE.UBOUND)RETURN
      IF(UOLD.LE.UBOUND)GO TO 20
      IF(USTEP.GT.STEPM)T=STEPM/USTEP
      RETURN
   20 T=(UBOUND+DELTA-UOLD)/USTEP
      RETURN
  100 CONTINUE
      RETURN
C
      END
