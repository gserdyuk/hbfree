c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




       SUBROUTINE INDSV(OM,P1,L1,P2,L2,P3,L3)
C
C     MODEL OF THREE MUTUALLY INDUCTIVELY COUPLED
C     TWO-PORT NETWORKS (LINEAR)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION OM,P1,P2,P3
      DOUBLE PRECISION R1, R2, R3
      DOUBLE PRECISION K1, K2, K3
      DOUBLE PRECISION LL1,LL2,LL3
      DOUBLE PRECISION M1, M2, M3
      DOUBLE PRECISION FLAG1,FLAG2,FLAG3
      DOUBLE PRECISION ZA1,ZA2,ZA3,ZB1,ZB2,ZB3
      DIMENSION P1(L1),P2(L2),P3(L3)
      COMMON/SUBC/YL(15,15),VJ(15)
      DOUBLE COMPLEX YL,VJ
      DOUBLE COMPLEX A,B,C,D,E,F
      DOUBLE COMPLEX Z1,Z2,Z3,ZM1,ZM2,ZM3,Z

C     PARAMETERS COMMON TO THE SCHEME: G1, G2, G3
C     INDIVIDUAL PARAMETERS:
C                Z1 = R1 + j*OM*LL1 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z1.
C                Z2 = R2 + j*OM*LL2 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z2.
C                Z3 = R3 + j*OM*LL3 - COMPLEX TOTAL
C     IMPEDANCE OF THE BRANCH CONTAINING TWO-PORT Z3.
C                M1 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z1 AND Z2.
C                M2 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z2 AND Z3.
C                M3 - MUTUAL INDUCTANCE OF TWO-PORTS
C                    Z3 AND Z1.
C                K1 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z1 AND Z2.
C                K2 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z2 AND Z3.
C                K3 - COUPLING COEFFICIENT OF TWO-PORTS
C                    Z3 AND Z1.
C            FLAG1 = 0. - TWO-PORTS Z1 AND Z2 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z1 AND Z2 CONNECTED OUT OF PHASE
C            FLAG2 = 0. - TWO-PORTS Z2 AND Z3 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z2 AND Z3 CONNECTED OUT OF PHASE
C            FLAG3 = 0. - TWO-PORTS Z3 AND Z1 CONNECTED IN PHASE
C                    1. - TWO-PORTS Z3 AND Z1 CONNECTED OUT OF PHASE
C

      G1=P2(1)
      G2=P2(1)
      G3=P2(1)
      IF(P3(1).EQ.0.0D0) P3(1)=1.D-11
      R1=P3(1)
      LL1=P3(2)
      IF(P3(3).EQ.0.0D0) P3(3)=1.D-11
      R2=P3(3)
      LL2=P3(4)
      IF(P3(5).EQ.0.0D0) P3(5)=1.D-11
      R3=P3(5)
      LL3=P3(6)
      K1=P3(7)
      K2=P3(8)
      K3=P3(9)
      FLAG1=P3(10)
      FLAG2=P3(11)
      FLAG3=P3(12)

      M1=K1*DSQRT(LL1*LL2)
      M2=K2*DSQRT(LL2*LL3)
      M3=K3*DSQRT(LL3*LL1)

      ZA1=OM*LL1
      Z1=DCMPLX(R1,ZA1)
      ZA2=OM*LL2
      Z2=DCMPLX(R2,ZA2)
      ZA3=OM*LL3
      Z3=DCMPLX(R3,ZA3)

      ZB1=OM*M1
      ZB2=OM*M2
      ZB3=OM*M3

      IF(FLAG1.EQ.0.0D0) ZM1=DCMPLX(0.0D0,ZB1)
      IF(FLAG1.EQ.1.0D0) ZM1=-DCMPLX(0.0D0,ZB1)
      IF(FLAG1.NE.0.0D0.AND.FLAG1.NE.1.0D0) GOTO 100

      IF(FLAG2.EQ.0.0D0) ZM2=DCMPLX(0.0D0,ZB2)
      IF(FLAG2.EQ.1.0D0) ZM2=-DCMPLX(0.0D0,ZB2)
      IF(FLAG2.NE.0.0D0.AND.FLAG2.NE.1.0D0) GOTO 100

      IF(FLAG3.EQ.0.0D0) ZM3=DCMPLX(0.0D0,ZB3)
      IF(FLAG3.EQ.1.0D0) ZM3=-DCMPLX(0.0D0,ZB3)
      IF(FLAG3.NE.0.0D0.AND.FLAG3.NE.1.0D0) GOTO 100

      WRITE(6,10) Z1,Z2,Z3,ZM1,ZM2,ZM3
   10 FORMAT(2X,'Z1 = ',E12.5,' j',E12.5/2X,'Z2 = ',E12.5,' j',E12.5/
     +    2X,'Z3 = ',E12.5,' j',E12.5/2X,'ZM1 = ',E12.5,' j',E12.5 /
     +   2X,'ZM2 = ',E12.5,' j',E12.5/2X,'ZM3 = ',E12.5,' j',E12.5)

      Z=Z1*Z2*Z3+2*ZM1*ZM2*ZM3-Z2*ZM3**2-Z3*ZM1**2-Z1*ZM2**2

      A=(Z2*Z3-ZM2**2)/Z
      B=(-ZM1*Z3+ZM2*ZM3)/Z
      C=(ZM1*ZM2-ZM3*Z2)/Z
      D=(Z1*Z3-ZM3**2)/Z
      E=(-ZM2*Z1+ZM1*ZM3)/Z
      F=(Z1*Z2-ZM1**2)/Z

      WRITE(6,15) Z,A,B,C,D,E,F
   15 FORMAT(2X,'Z = ',E12.5,' j',E12.5/2X,'A = ',E12.5,' j',E12.5/
     +  2X,'B = ',E12.5,' j',E12.5/2X,'C = ',E12.5,' j',E12.5/       2X,
     +'D = ',E12.5,' j',E12.5/2X,'E = ',E12.5,' j',E12.5/       2X,'F =
     +',E12.5,' j',E12.5)

      YL(1,1)=A+G1
      YL(1,2)=-A
      YL(1,3)=B-G1
      YL(1,4)=-B
      YL(1,5)=C
      YL(1,6)=-C
      YL(2,1)=-A
      YL(2,2)=A+G3
      YL(2,3)=-B
      YL(2,4)=B
      YL(2,5)=-C
      YL(2,6)=C-G3
      YL(3,1)=B-G1
      YL(3,2)=-B
      YL(3,3)=D+G1
      YL(3,4)=-D
      YL(3,5)=E
      YL(3,6)=-E
      YL(4,1)=-B
      YL(4,2)=B
      YL(4,3)=-D
      YL(4,4)=D+G2
      YL(4,5)=-E-G2
      YL(4,6)=E
      YL(5,1)=C
      YL(5,2)=-C
      YL(5,3)=E
      YL(5,4)=-E-G2
      YL(5,5)=F+G2
      YL(5,6)=-F
      YL(6,1)=-C
      YL(6,2)=C-G3
      YL(6,3)=-E
      YL(6,4)=E
      YL(6,5)=-F
      YL(6,6)=F+G3

      WRITE(6,20) ((N1,N2,YL(N1,N2),N2=1,6),N1=1,6)
  20  FORMAT(2X,'INDSV',2X,'YL(',I2,',',I2,')=',E12.5,2X,E12.5)

      RETURN

100   CONTINUE
      WRITE (6,30)
30    FORMAT(2X,52('$')/      2X,'$$                E R R O R
     +             $$'/       2X,'$$    IN THE INPUT DATA. 10-12 INDIVID
     +UAL          $$'/      2X,'$$    PARAMETERS FOR THE INDUCTIVELY-ÛO
     +UPLED       $$'/      2X,'$$    TWO-PORT NETWORKS CAN ONLY HAVE
     +           $$'/      2X,'$$    THE FOLLOWING VALUES:
     +          $$'/      2X,'$$    0. - TWO-PORTS CONNECTED IN PHASE.
     +         $$'/      2X,'$$    1. - TWO-PORTS CONNECTED OUT OF PHASE
     +.       $$'/      2X,52('$'))

      RETURN

      END
