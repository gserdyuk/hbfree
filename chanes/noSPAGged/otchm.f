c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE OTCHM(Y,KOL)
c$LARGE: Y
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION Y(1)
      KY2=KOL+KOL
      DO 5 I=2,KY2,2
    5 Y(I)=0.0D0
C      PRINT 10, KY2,(Y(J),J=1,KY2)
C   10 FORMAT (2x,'OTCHM: KY2=',I4/2X,'Y=',3X,6(E13.6))
      RETURN
C     DEBUG SUBTRACE
      END
