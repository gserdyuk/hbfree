!*==OTCHM.for processed by SPAG 8.04DB 09:58  3 May 2025
!!SPAG Open source Personal, Educational or Academic User Hobbyist  NON-COMMERCIAL USE - Not for use on proprietary or closed source code
c
c Copyright (c) 1996-2004 by Gennady Serdyuk.  All rights reserved.
c gserdyuk@mail.ru
c
c Released under GPL v 2.0
c




      SUBROUTINE otchm(Y,Kol)
c$LARGE: Y
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      INTEGER i , Kol , ky2
      DOUBLE PRECISION Y
      DIMENSION Y(1)
      ky2 = Kol + Kol
      DO i = 2 , ky2 , 2
         Y(i) = 0.0D0
      ENDDO
C      PRINT 10, KY2,(Y(J),J=1,KY2)
C   10 FORMAT (2x,'OTCHM: KY2=',I4/2X,'Y=',3X,6(E13.6))
C     DEBUG SUBTRACE
      END
