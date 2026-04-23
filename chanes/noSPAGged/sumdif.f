c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE SUMDIF(Q,IR,IS,SUM,DIFF)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON /BLW1/   W, W1
      DOUBLE PRECISION            W(20),W1(200)
      COMMON /BLW2/   WR, WS
      INTEGER         WR(20,20), WS(20,20)
      DOUBLE COMPLEX         Q(1), SUM, DIFF, QIZS
C
C
      IZS=WS(IR,IS)
      QIZS=Q(IZS)
      IZR=WR(IR,IS)
      IF(IZR.LT.0) GO TO 10
      SUM=Q(IZR)
      GO TO 20
   10 IZR=-IZR
      SUM=DCONJG(Q(IZR))
   20 DIFF=SUM-QIZS
      IF(W(IS).NE.0.D0) SUM=SUM+QIZS
      RETURN
C     DEBUG INIT(Q,IR,IS,SUM,DIFF),SUBTRACE
      END
