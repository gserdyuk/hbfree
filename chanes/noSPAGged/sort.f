c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE SORT (MN,KN)
C
C  *****  SORTING OF THE MN ARRAY
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MN(2,KN)
C
      IF (KN.LT.2) RETURN
      K1=KN-1
      DO 15 I=1,K1
      M=MN(1,I)
      N=MN(2,I)
      KRIT=1024*M+N
      IF (N .LT.0) KRIT =KRIT +256
      J1=I
      I1=I+1
      DO 10 J=I1,KN
      M1=MN(1,J)
      N1=MN(2,J)
      KRIT1=1024*M1+N1
      IF (N1.LT.0) KRIT1=KRIT1+256
      IF (KRIT.LE.KRIT1) GO TO 10
      KRIT=KRIT1
      J1=J
   10 CONTINUE
      IF (J1.EQ.I) GO TO 15
      MNM=MN(1,I)
      MN(1,I)=MN(1,J1)
      MN(1,J1)=MNM
      MNM=MN(2,I)
      MN(2,I)=MN(2,J1)
      MN(2,J1)=MNM
   15 CONTINUE
      RETURN
      END
