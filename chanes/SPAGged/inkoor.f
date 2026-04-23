c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE INKOOR (MAXSYN,MNMAX,F1,F2,*,*)
C
C  *****  ORDERING OF ARRAYS MN AND MN1,
C  *****  FILLING THE COORDINATE ARRAYS,
C  *****  OBTAINING THE EXPANDED FREQUENCY GRID,
C  *****  FILLING THE ARRAYS WR AND WS.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'funcsize.i'
      COMMON/BLW1/     W,W1
      DOUBLE PRECISION             W(20),W1(200)
      COMMON /BLW2/    WR,WS
      INTEGER          WR(20,20),WS(20,20)
      COMMON /BLK1/    KR,KR1,KC,KC1,NNR,NNR1,MN,MN1
      INTEGER          KR(20),KC(10),NNR(10)
      INTEGER          KR1(200),KC1(20),NNR1(20)
      INTEGER          MN(2,20), MN1(2,200)
      COMMON /BLK2/    KNC,KNR,KN,KNR1,KN1
      COMMON / MEP /   MEPHF, FLGMNW
      INTEGER          FLGMNW
      INTEGER          B1DIM,B2DIM
      DATA             IR/10/
C
C take sizes from 'funcsize.i'
      B1DIM = B1_SIZE
      B2DIM = B2_SIZE
      NDIM  = MAXKN
      NDIM1 = MAXKN1

      PI=4*DATAN(1.0D0)
      OMEGA1=2.0D0*PI*F1
      OMEGA2=2.0D0*PI*F2
      DO 10 I=1,KN
      IF (MN(1,I).GT.0) GO TO 10
      IF (MN(1,I).EQ.0.AND.MN(2,I).GE.0) GO TO 10
      MN(1,I)=-MN(1,I)
      MN(2,I)=-MN(2,I)
   10 CONTINUE
      FLGMNW=0
      MEPHF=2
      IF(KN.EQ.1) GO TO 1065
      CALL SORT (MN,KN)




C  *****  CHANGE OF FREQUENCIES F1 AND F2, IF F2 = 0

      DO 1010 I=1,KN
      IF (MN(2,I).NE.0) GO TO 1030
 1010 CONTINUE
      DO 1020 I=1,KN
      MNM=MN(1,I)
      MN(1,I)=MN(2,I)
 1020 MN(2,I)=MNM
      OMEGA =OMEGA1
      OMEGA1=OMEGA2
      OMEGA2=OMEGA
      MEPHF=1
      FLGMNW=1
C  *****  EXCLUSION OF IDENTICAL ROWS OF THE MATRIX
C  *****  AND COMPRESSION OF THE LIST MN AFTER EXCLUSION
C
 1030 I=2
 1040 IF(I.GT.KN) GO TO 1065
      IF    (MN(1,I).EQ.MN(1,I-1) .AND. MN(2,I).EQ.MN(2,I-1)) GO TO 1050
      I=I+1
      GO TO 1040
 1050 KN=KN-1
      IF(I.GT.KN) GO TO 1065
      DO 1060 J=I,KN
      MN(1,J)=MN(1,J+1)
 1060 MN(2,J)=MN(2,J+1)
      GO TO 1040
C  *****  ADDITION OF THE CONSTANT COMPONENT,
C  *****  IF IT IS ABSENT IN THE LIST MN
 1065 IF(MN(1,1).EQ.0.AND.MN(2,1).EQ.0) GO TO 1100
      DO 1080 I=1,KN
      MN(1,KN-I+2)=MN(1,KN-I+1)
 1080 MN(2,KN-I+2)=MN(2,KN-I+1)
      MN(1,1)=0
      MN(2,1)=0
      KN=KN+1
 1100 IF (KN.EQ.1) MEPHF=0
      DO 15 I=1,KN
   15 W (I)=MN (1,I)*OMEGA1+MN (2,I)*OMEGA2
      KN1=0
      DO 70 I=1,KN
      I1=MN(1,I)
      I2=MN(2,I)
      DO 70 J=1,KN
      J1=MN(1,J)
      J2=MN(2,J)
      MP=1
      DO 60 L=1,2
      IPQ1=I1+MP*J1
      IPQ2=I2+MP*J2
      IF (IPQ1) 30,35,45
   30 IPQ1=-IPQ1
      GO TO 40
   35 IF(IPQ2.GE.0) GO TO 45
   40 IPQ2=-IPQ2
   45 IF (KN1.EQ.0) GO TO 55
      DO 50 K=1,KN1
      IF (IPQ1.EQ.MN1(1,K).AND.IPQ2.EQ.MN1(2,K)) GO TO 60
   50 CONTINUE
   55 KN1=KN1+1
      IF (KN1.GT.NDIM1) GO TO 200
      MN1(1,KN1)=IPQ1
      MN1(2,KN1)=IPQ2
   60 MP=-1
   70 CONTINUE
      CALL SORT (MN1,KN1)
      DO 73 I=1,KN1
   73 W1(I)=MN1(1,I)*OMEGA1+MN1(2,I)*OMEGA2
C  *****  CALCULATION OF THE NUMBER OF NONZERO ROWS
C  *****  ON THE EXPANDED FREQUENCY GRID - KNR1
      KNR1=1
      IF (KN.LT.2) GO TO 2010
      DO 2000 I=2,KN1
      IF (MN1(1,I).EQ.MN1(1,I-1)) GO TO 2000
      KNR1=KNR1+1
 2000 CONTINUE
C  *****  SELECTION OF LIMITS MIN AND MAX FOR KNC
 2010 MNMAX=0
      ALOG2=DLOG10(2.D0)
      DO 2020 I=1,KN
 2020 MNMAX=MAX0(MNMAX,MN(1,I),IABS(MN(2,I)))
      KNCMIN=MAX0(MNMAX*4,16)
      KNCMIN=2**(IDINT(DLOG10(KNCMIN-0.5D0)/ALOG2)+1)
      KNCMAX=MIN0(B1DIM/MAXSYN,B2DIM/(MAXSYN*KNR1),128)
      KNCMAX=2**(IDINT(DLOG10(KNCMAX+0.5D0)/ALOG2))
      IF (KNCMIN.GT.KNCMAX) GO TO 250
      IF (KNC.GT.KNCMIN.AND.KNC.LT.KNCMAX) KNC=2**(IDINT(DLOG10(KNC-0.5D
     +0)/ALOG2)+1)
      IF (KNC.LE.KNCMIN) KNC=KNCMIN
      IF (KNC.GE.KNCMAX) KNC=KNCMAX
C  *****  FILLING THE COORDINATE ARRAYS
      CALL KOORD (MN,KR,KC,NNR,KNR,KN,KNC,IR,NDIM)

C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
C                          (TOTAL 3 CHANGES)
      ivald= krdchk (MN,KR,KC,NNR,KNR,KN,KNC,IR,NDIM)
      if(ivald.ne.0) then
                     write(6,5010) ivald
                     print 5010,ivald
      endif
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

C  *****  FILLING THE COORDINATE ARRAYS
C  *****  FOR THE EXPANDED FREQUENCY GRID
      CALL KOORD (MN1,KR1,KC1,NNR1,KNR1,KN1,KNC,NDIM,NDIM1)
C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
      ivald1=krdchk (MN1,KR1,KC1,NNR1,KNR1,KN1,KNC,NDIM,NDIM1)
      if(ivald1.ne.0) then
                     write(6,5020)ivald1
                     print 5020,ivald1
      endif
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NNR(KNR).EQ.1.AND.MEPHF.EQ.2) MEPHF=1
C  *****  FILLING THE ARRAYS WR AND WS
      DO 120 I=1,KN
      I1=MN(1,I)
      I2=MN(2,I)
      DO 120 J=1,KN
      J1=MN(1,J)
      J2=MN(2,J)
      MP=1
      DO 110 L=1,2
      IPQ1=I1+MP*J1
      IPQ2=I2+MP*J2
      KK=1
      IF (IPQ1) 75,80,90
   75 IPQ1=-IPQ1
      GO TO 85
   80 IF (IPQ2.GE.0) GO TO 90
   85 IPQ2=-IPQ2
      KK=-1
   90 DO 95 K=1,KN1
      IF (IPQ1.EQ.MN1(1,K) .AND. IPQ2.EQ.MN1(2,K)) GO TO 100
   95 CONTINUE
  100 K=K*KK
      IF (L.EQ.1) WS(I,J)=K
      IF (L.EQ.2) WR(I,J)=K
  110 MP=-1
  120 CONTINUE
      RETURN
C     DEBUG INIT(W,OMEGA1,OMEGA2,MN)
C----------------------------------------------------------------------
C  *****  THE SPREAD OF MN1 EXCEEDS THE MAXIMUM
  200 CONTINUE
      RETURN 1
C----------------------------------------------------------------------
C  *****  THE SPREAD OF THE TRANSFORMATION OF THE FOURIER
C  *****  EXCEEDS THE UPPER LIMIT VALUE
C
  250 MNMAX=KNCMAX/4

      RETURN 2
C++++++++++++++++++++++++  CHANGES FROM 12.05.91  SERDYUK G.V.
 5010 format(2x,' INKOOR: ERROR OF REPRESENTATION OF MAIN GRID.',
     +   ' CODE =',i3)
 5020 format(2x,' INKOOR: ERROR OF REPRESENTATION OF AUX  GRID.',
     +   ' CODE =',i3)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      END
