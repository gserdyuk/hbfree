c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE KOORD(MN,KR,KC,NNR,KNR,KN,KNC,IR1,IR2)
C
C  *****  FILLING THE COORDINATE ARRAYS
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER    KR(IR2),KC(IR1),NNR(IR1)
      INTEGER    MN(2,IR2)
C
      DO 10 I=1,IR1
      KC(I)=0
   10 NNR(I)=0
      DO 15 I=1,IR2
   15 KR(I)=0
      KNR=0
      MNR=-1
      DO 40 I=1,KN
      N=MN(2,I)
      IF (N.LT.0) N=N+KNC
      KR(I)=N+1
      M=MN(1,I)
      IF (M.EQ.MNR) GO TO 17
      KNR=KNR+1
      NNR(KNR)=M+1
      MNR=M
   17 CONTINUE
      IF (I.GT.1) GO TO 20
      MS=1
      M1=M
      J=1
C ++++++++++++++++++++++++++++ CHANGE FROM 12.05.91
C                             TOTAL 1 PIECE.
      IF (I.EQ.KN) then
           GO TO 30
      else
                   GO TO 40
      endif
   20 continue
      IF (M.NE.M1.and.i.lt.kn) then
                   GO TO 30
      elseif (.not.M.NE.M1.and.i.lt.kn)then
                                       ms=ms+1
                                       goto 40
      elseif (M.NE.M1.and..not.i.lt.kn)then
                                       kc(j)=ms
                                       kc(j+1)=1
                                       goto 40
      else
                                       ms=ms+1
                                       goto 30
      endif
C +++++++++++++++++++++++++++++++++++++++++++++++
   30 KC(J)=MS
      MS=1
      M1=M
      J=J+1
   40 CONTINUE
      RETURN
      END
C ++++++++++++++++++++++++++++++++++++++++++++++++
C  A POSSIBLE VARIANT OF THE SUBROUTINE/PROGRAM:
C                       (IF THERE IS TIME - TRY IT)
C
C      SUBROUTINE KOORD(MN,KR,KC,NNR,KNR,KN,KNC,IR1,IR2)
CC
CC  *****  FILLING THE COORDINATE ARRAYS
CC
C      INTEGER    KR(IR2),KC(IR1),NNR(IR1)
C      INTEGER    MN(2,IR2)
C
C      DO 10 I=1,IR1
C      KC(I)=0
C   10 NNR(I)=0
C      DO 15 I=1,IR2
C   15 KR(I)=0
C
C      KNR=1
C      MN11=MN11
C      IF(MN11.LT.0)MN11=MN11+KNC
C      NNR(1)=MN11+1
C
C      DO 40 I=1,KN
C
C      MN1I=MN(1,I)
C      MN2I=MN(2,I)
C      IF(MN1I.LT.0)MN1I=MN1I+KNC
C      IF(MN2I.LT.0)MN2I=MN2I+KNC
C
C      IF( NNR(KNR).NE.(MN1I+1)) THEN
C                                KNR=KNR+1
C                                NNR(KNR)=MN1I+1
C      ENDIF
C
C      KC(KNR)=KC(KNR)+1
C      KR(I)=MN2I+1
C 40   CONTINUE
C
C      RETURN
C      END
C
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
C                             Serdyuk G.V.   12.05.91.
C +++++++++++++++++++++++++++++++++++++++++++++++++++++++
